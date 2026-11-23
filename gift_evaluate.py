"""gift_evaluate.py -- evaluation of GIFT-64 with its original key schedule and with GAKS-GIFT.
Validation keys (seed 777) decide the warm-up count as for DES; evaluation keys use seed 20261122.
Writes results/gift/eval.json."""
import json, math, sys
import numpy as np
from scipy.stats import norm
from gaks_core import U, popcount
from gift_core import *
from ga_design import bits_of

NAME = sys.argv[1] if len(sys.argv) > 1 else "nb2_s7"
g = json.load(open(f"results/gift/ga_gift_{NAME}.json"))
PRM = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"])
CELLS = 128 * ROUNDS * 32
ZB = float(norm.isf(0.01 / 2 / CELLS))
KS = {"GIFT-64": lambda h, l: ks_gift(h, l), "GAKS-GIFT": lambda h, l: ks_gaks_gift(h, l, PRM)}

def keys(rng, n):
    return rng.integers(0, 2**63, n, dtype=np.uint64) * U(2) + rng.integers(0, 2, n, dtype=np.uint64), \
           rng.integers(0, 2**63, n, dtype=np.uint64) * U(2) + rng.integers(0, 2, n, dtype=np.uint64)

def flipped(H, L):
    FH = np.repeat(H[:, None], 128, 1).copy(); FL = np.repeat(L[:, None], 128, 1).copy()
    for b in range(64):
        FL[:, b] ^= U(1) << U(b); FH[:, 64 + b] ^= U(1) << U(b)
    return FH, FL

def diffusion(ksf, H, L):
    n = len(H); FH, FL = flipped(H, L)
    base = ksf(H, L); x = ksf(FH.ravel(), FL.ravel()).reshape(n, 128, ROUNDS) ^ base[:, None, :]
    hd = popcount(x) / 32; p = bits_of(x, 32).mean(0); d = np.abs(p - 0.5)
    sig = (d * 2 * np.sqrt(n) > ZB).sum(axis=(0, 2))
    return dict(per_round=hd.mean(axis=(0, 1)).tolist(), mean=float(hd.mean()), sac=float(d.mean()),
                sig_per_round=sig.tolist(), det_cells=int((d == 0.5).sum())), x, base

R = {"ga_params": PRM.to_dict(), "bonferroni_z": ZB}
# validation rule
VK = keys(np.random.default_rng(777), 500)
R["validation"] = []
while True:                       # same rule as for DES, applied repeatedly: add warm-up steps until clean
    vd, _, _ = diffusion(lambda h, l: ks_gaks_gift(h, l, PRM), *VK)
    R["validation"].append({"w": PRM.w, "sig_per_round": vd["sig_per_round"]})
    print("validation w", PRM.w, vd["sig_per_round"][:4], flush=True)
    if (sum(vd["sig_per_round"][1:]) == 0 and vd["sig_per_round"][0] <= 5) or PRM.w >= 12:
        break
    PRM.w += 1
R["final_w"] = PRM.w; R["sboxes"] = PRM.sbox_count(); R["params"] = PRM.to_dict()
rng = np.random.default_rng(20261122)
H, L = keys(rng, 500)
R["ideal_sac_ref"] = float(np.abs(rng.integers(0, 2, (500, 128, ROUNDS, 32)).mean(0) - 0.5).mean())
for name, f in KS.items():
    d, x, base = diffusion(f, H, L)
    # E4: max probability of the most frequent round-key difference over single-bit key differences
    worst = np.zeros(ROUNDS)
    for b in range(128):
        for r in range(ROUNDS):
            _, c = np.unique(x[:, b, r], return_counts=True); worst[r] = max(worst[r], c.max() / len(H))
    d["rk_max_prob_per_round"] = worst.tolist()
    R[name] = d
    print(name, round(d["mean"], 4), round(d["sac"], 4), d["sig_per_round"][:4], d["det_cells"], worst.max())
# cipher key avalanche (round-reduced)
H3, L3 = H[:200], L[:200]; FH, FL = flipped(H3, L3); pt = rng.integers(0, 2**63, 200, dtype=np.uint64)
rnd = [GAKSGiftParams.from_bits(np.random.default_rng(600 + i).integers(0, 2, GAKSGiftParams.nbits(PRM.nblocks)), PRM.nblocks) for i in range(10)]
def cav(f):
    _, r0 = gift_encrypt(pt, f(H3, L3), record=True); _, r1 = gift_encrypt(np.repeat(pt, 128), f(FH.ravel(), FL.ravel()), record=True)
    return (popcount(r1 ^ np.repeat(r0, 128, axis=0)) / 64).mean(axis=0)
for name, f in KS.items():
    R[name]["cipher_av"] = cav(f).tolist()
curves = np.array([cav(lambda h, l, q=q: ks_gaks_gift(h, l, q)) for q in rnd])
R["random_designs_cipher_av"] = dict(mean=curves.mean(0).tolist(), min=curves.min(0).tolist(), max=curves.max(0).tolist())
for name in KS:
    R[name]["rounds_to_049"] = int(next(i + 1 for i, v in enumerate(R[name]["cipher_av"]) if v >= 0.49))
# BLR linearity
a, b, c = keys(rng, 3000), keys(rng, 3000), keys(rng, 3000)
for name, f in KS.items():
    lhs = f(a[0] ^ b[0] ^ c[0], a[1] ^ b[1] ^ c[1]); rhs = f(*a) ^ f(*b) ^ f(*c)
    R[name]["blr_full"] = float(np.mean(np.all(lhs == rhs, axis=1)))
    R[name]["blr_bits"] = float((1 - popcount(lhs ^ rhs).sum(axis=1) / (32 * ROUNDS)).mean())
# round-key bit dependence (K1..K8)
hh, ll = keys(rng, 20000)
def bitcorr(rk):
    B = bits_of(rk, 32).reshape(rk.shape[0], -1).astype(np.float32); B -= B.mean(0); B /= B.std(0) + 1e-9
    C = (B.T @ B) / B.shape[0]; blk = np.kron(np.eye(rk.shape[1]), np.ones((32, 32))).astype(bool); C[blk] = 0
    return float(np.abs(C).max())
for name, f in KS.items():
    R[name]["interkey_bitcorr"] = bitcorr(f(hh, ll)[:, :8])
R["interkey_noise"] = bitcorr(rng.integers(0, 2**32, (20000, 8), dtype=np.uint64))
# period of the original (linear, bit-permutation) key-state update
perm = []
for bpos in range(128):
    h = as_u64((1 << (bpos - 64)) if bpos >= 64 else 0); l = as_u64((1 << bpos) if bpos < 64 else 0)
    h, l = gift_ks_update(h, l); vv = (int(h[0]) << 64) | int(l[0]); perm.append(vv.bit_length() - 1)
seen, cyc = set(), []
for i in range(128):
    if i in seen: continue
    j, n = i, 0
    while j not in seen: seen.add(j); j = perm[j]; n += 1
    cyc.append(n)
R["GIFT-64"]["ks_period"] = math.lcm(*cyc)
json.dump(R, open("results/gift/eval.json", "w"), indent=1, default=float)
print("final_w", PRM.w, "sboxes", PRM.sbox_count(), "period", R["GIFT-64"]["ks_period"])
print({k: (R[k]["rounds_to_049"], R[k]["blr_bits"], R[k]["interkey_bitcorr"]) for k in KS}, R["interkey_noise"])
