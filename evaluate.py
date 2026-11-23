"""
evaluate.py -- Security/diffusion evaluation of standard DES / PRESENT-80 key schedules versus
the GA-designed schedules (GAKS). Uses evaluation keys drawn from seeds disjoint from the design set.
Writes results/eval.json.
"""
import json, time
import numpy as np
from gaks_core import *
from ga_design import bits_of

rng = np.random.default_rng(20261122)   # fresh evaluation seed (777 was used as validation set)
gd = json.load(open("results/ga_des.json"))
gp = json.load(open("results/ga_present.json"))
PRM = final_des_params()
PRM_GA = GAKSParams.from_bits(gd["chromosome"])   # raw GA output (w=6)
PRMP = GAKSPresentParams.from_bits(gp["chromosome"])
RES = {"params_des": PRM.to_dict(), "params_present": PRMP.to_dict()}

def rand_keys(n):
    return rng.integers(0, 2**63, n, dtype=np.uint64) * U(2) + rng.integers(0, 2, n, dtype=np.uint64)

KS = {"DES": ks_des, "GAKS-DES": lambda k: ks_gaks(k, PRM)}
rand_designs = [GAKSParams.from_bits(np.random.default_rng(100 + i).integers(0, 2, GAKSParams.NBITS)) for i in range(10)]
FLIPS = np.array([U(1) << U(b) for b in EFFECTIVE_BITS], dtype=np.uint64)

# ------------------------------------------------ E1 correctness & determinism
t = time.time()
keys = rand_keys(10000); pts = rand_keys(10000)
sk = KS["GAKS-DES"](keys)
assert np.array_equal(sk, KS["GAKS-DES"](keys)), "non-deterministic"
ct = des_encrypt(pts, sk)
assert np.array_equal(des_decrypt(ct, sk), pts)
RES["E1"] = {"pairs": 10000, "roundtrip_ok": True, "deterministic": True}
print("E1 ok", time.time() - t)

# ------------------------------------------------ E2/E3 key-schedule avalanche + SAC
def ks_diffusion(ksf, N, flips, nbits, R, keygen):
    base, fl = keygen(N, ksf, flips)
    x = fl ^ base[:, None, :]
    hd = popcount(x) / nbits                        # (N,F,R)
    p = bits_of(x, nbits).mean(axis=0)              # (F,R,nbits)
    return hd, p

def des_keygen(N, ksf, flips):
    k = rand_keys(N)
    base = ksf(k)
    fl = ksf((k[:, None] ^ flips[None, :]).ravel()).reshape(N, len(flips), 16)
    return base, fl

N2 = 1000
from scipy.stats import norm as _norm
ZB = float(_norm.isf(0.01 / 2 / (56 * 16 * 48)))     # two-sided, Bonferroni over all 43,008 SAC cells
E2 = {"bonferroni_z": ZB}
for name, ksf in list(KS.items()):
    hd, p = ks_diffusion(ksf, N2, FLIPS, 48, 16, des_keygen)
    E2[name] = {"per_round_mean": hd.mean(axis=(0, 1)).tolist(), "per_round_std": hd.std(axis=(0, 1)).tolist(),
                "overall_mean": float(hd.mean()), "seq768_mean": float(hd.mean(axis=2).mean()),
                "sac_mean_abs_dev": float(np.abs(p - 0.5).mean()), "sac_max_abs_dev": float(np.abs(p - 0.5).max()),
                "sac_frac_within_0.05": float((np.abs(p - 0.5) <= 0.05).mean()),
                "sac_per_round_mean": np.abs(p - 0.5).mean(axis=(0, 2)).tolist(),
                "sac_per_round_max": np.abs(p - 0.5).max(axis=(0, 2)).tolist(),
                "sac_per_round_frac_gt_0.05": (np.abs(p - 0.5) > 0.05).mean(axis=(0, 2)).tolist(),
                "deterministic_cells": int((np.abs(p - 0.5) == 0.5).sum()),
                "significant_cells_per_round_bonferroni_1pct": (np.abs(p - 0.5) * 2 * np.sqrt(N2) > ZB).sum(axis=(0, 2)).tolist(),
                "cells_per_round": int(p.shape[0] * p.shape[2])}
    print("E2", name, E2[name]["overall_mean"], E2[name]["sac_mean_abs_dev"])
# reference: an ideal random function with the same sample size
ref = rng.integers(0, 2, (N2, 56, 16, 48)).mean(axis=0)
E2["ideal_random_reference"] = {"sac_mean_abs_dev": float(np.abs(ref - 0.5).mean()), "sac_max_abs_dev": float(np.abs(ref - 0.5).max()),
                                "sac_frac_within_0.05": float((np.abs(ref - 0.5) <= 0.05).mean())}
rd = []
for prm in rand_designs:
    hd, p = ks_diffusion(lambda k: ks_gaks(k, prm), 300, FLIPS, 48, 16, des_keygen)
    rd.append({"min_round_mean": float(hd.mean(axis=(0, 1)).min()), "sac_mean_abs_dev": float(np.abs(p - 0.5).mean()),
               "sboxes": prm.sbox_count(), "per_round_mean": hd.mean(axis=(0, 1)).tolist()})
E2["random_designs_N300"] = rd
RES["E2"] = E2

# ------------------------------------------------ E4 related-key predictability
N4 = 1000
E4 = {}
for name, ksf in KS.items():
    k = rand_keys(N4)
    base = ksf(k)
    worst_round, worst_seq = 0.0, 0.0
    per_round_worst = np.zeros(16)
    for f in FLIPS:
        d = ksf(k ^ f) ^ base                               # (N,16)
        for r in range(16):
            _, c = np.unique(d[:, r], return_counts=True)
            pr = c.max() / N4
            per_round_worst[r] = max(per_round_worst[r], pr)
        _, c = np.unique(d, axis=0, return_counts=True)
        worst_seq = max(worst_seq, c.max() / N4)
    E4[name] = {"max_prob_subkey_diff_per_round": per_round_worst.tolist(),
                "max_prob_over_all_rounds": float(per_round_worst.max()), "max_prob_full_sequence": float(worst_seq)}
    print("E4", name, per_round_worst.max(), worst_seq)
RES["E4"] = E4

# ------------------------------------------------ E5 complementation property
N5 = 10000
k = rand_keys(N5); p = rand_keys(N5)
E5 = {}
for name, ksf in KS.items():
    c1 = des_encrypt(p, ksf(k)); c2 = des_encrypt(~p, ksf(~k))
    E5[name] = {"holds_fraction": float(np.mean(c2 == ~c1)), "trials": N5}
print("E5", E5)
RES["E5"] = E5

# ------------------------------------------------ E6 weak / semi-weak keys
WEAK = [0x0101010101010101, 0xFEFEFEFEFEFEFEFE, 0xE0E0E0E0F1F1F1F1, 0x1F1F1F1F0E0E0E0E]
SEMI = [(0x01FE01FE01FE01FE, 0xFE01FE01FE01FE01), (0x1FE01FE00EF10EF1, 0xE01FE01FF10EF10E),
        (0x01E001E001F101F1, 0xE001E001F101F101), (0x1FFE1FFE0EFE0EFE, 0xFE1FFE1FFE0EFE0E),
        (0x011F011F010E010E, 0x1F011F010E010E01), (0xE0FEE0FEF1FEF1FE, 0xFEE0FEE0FEF1FEF1)]
p = rand_keys(1000)
E6 = {}
for name, ksf in KS.items():
    w_rows = []
    for kw in WEAK:
        s = ksf(as_u64(kw))
        s2 = np.repeat(s, 1000, axis=0)
        inv = float(np.mean(des_encrypt(des_encrypt(p, s2), s2) == p))
        w_rows.append({"key": f"{kw:016X}", "distinct_subkeys": int(len(set(s[0].tolist()))), "E_k(E_k(p))=p_fraction": inv})
    s_rows = []
    for a, b in SEMI:
        sa = np.repeat(ksf(as_u64(a)), 1000, axis=0); sb = np.repeat(ksf(as_u64(b)), 1000, axis=0)
        frac = float(np.mean(des_encrypt(des_encrypt(p, sa), sb) == p))
        s_rows.append({"pair": f"{a:016X}/{b:016X}", "distinct_subkeys_a": int(len(set(ksf(as_u64(a))[0].tolist()))),
                       "E_k2(E_k1(p))=p_fraction": frac})
    E6[name] = {"weak": w_rows, "semi_weak": s_rows}
# random-key check: any repeated subkey within a schedule?
k = rand_keys(100000)
for name, ksf in KS.items():
    s = ksf(k)
    srt = np.sort(s, axis=1)
    E6[name]["random_keys_with_repeated_subkey"] = int(np.any(srt[:, 1:] == srt[:, :-1], axis=1).sum())
    E6[name]["random_keys_tested"] = 100000
k = rand_keys(100000)
for name, ksf in KS.items():
    s = ksf(k)
    E6[name]["random_keys_all_subkeys_equal"] = int(np.all(s == s[:, :1], axis=1).sum())
    E6[name]["random_keys_involutive_schedule"] = int(np.all(s == s[:, ::-1], axis=1).sum())
print("E6 done")
RES["E6"] = E6

# ------------------------------------------------ E7 round-reduced key avalanche (full cipher)
N7 = 300
E7 = {}
k = rand_keys(N7); p = rand_keys(N7)
kk = (k[:, None] ^ FLIPS[None, :]).ravel(); pp = np.repeat(p, 56)
def cipher_av(ksf):
    _, r0 = des_encrypt(p, ksf(k), record=True)
    _, r1 = des_encrypt(pp, ksf(kk), record=True)
    return popcount(r1 ^ np.repeat(r0, 56, axis=0)) / 64       # (N*56,16)
for name, ksf in KS.items():
    hd = cipher_av(ksf)
    E7[name] = {"per_round_mean": hd.mean(axis=0).tolist(), "per_round_std": hd.std(axis=0).tolist()}
    print("E7", name, np.round(hd.mean(axis=0)[:6], 3))
rnd20 = [GAKSParams.from_bits(np.random.default_rng(300 + i).integers(0, 2, GAKSParams.NBITS)) for i in range(20)]
curves = np.array([cipher_av(lambda x, q=q: ks_gaks(x, q)).mean(axis=0) for q in rnd20])
E7["random-designs-20"] = {"per_round_mean": curves.mean(axis=0).tolist(), "per_round_min": curves.min(axis=0).tolist(),
                           "per_round_max": curves.max(axis=0).tolist()}
RES["E7"] = E7

# ------------------------------------------------ E8 BLR linearity test on the key schedule
N8 = 10000
a, b, c = rand_keys(N8), rand_keys(N8), rand_keys(N8)
E8 = {}
for name, ksf in KS.items():
    lhs = ksf(a ^ b ^ c); rhs = ksf(a) ^ ksf(b) ^ ksf(c)
    eq_bits = 1 - popcount(lhs ^ rhs).sum(axis=1) / 768
    E8[name] = {"full_sequence_linear_fraction": float(np.mean(np.all(lhs == rhs, axis=1))),
                "mean_fraction_of_bits_satisfying": float(eq_bits.mean())}
print("E8", E8)
RES["E8"] = E8

# ------------------------------------------------ E10 cost
def timeit(f, n=20000, rep=5):
    k = rand_keys(n); best = 1e9
    for _ in range(rep):
        t = time.perf_counter(); f(k); best = min(best, time.perf_counter() - t)
    return n / best
RES["E10"] = {"DES_keys_per_s_vectorised": timeit(ks_des), "GAKS_keys_per_s_vectorised": timeit(KS["GAKS-DES"]),
              "DES_sboxes": 0, "GAKS_sboxes": PRM.sbox_count(), "GAKS_steps": PRM.w + 16}
print("E10", RES["E10"])

# ------------------------------------------------ E11 PRESENT-80 generalisation
def rand40(n): return rng.integers(0, 2**40, n, dtype=np.uint64)
KSP = {"PRESENT-80": lambda h, l: ks_present(h, l), "GAKS-PRESENT": lambda h, l: ks_gaks_present(h, l, PRMP)}
def present_flips(h, l):
    H, Lo = [], []
    for bpos in range(80):
        if bpos >= 40: H.append(h ^ U(1 << (bpos - 40))); Lo.append(l)
        else: H.append(h); Lo.append(l ^ U(1 << bpos))
    return np.stack(H, 1), np.stack(Lo, 1)
N11 = 500
h, l = rand40(N11), rand40(N11)
fh, fl_ = present_flips(h, l)
E11 = {}
pt = rand_keys(N11)
for name, ksf in KSP.items():
    base = ksf(h, l)
    fl = ksf(fh.ravel(), fl_.ravel()).reshape(N11, 80, 32)
    x = fl ^ base[:, None, :]
    hd = popcount(x) / 64
    pm = bits_of(x, 64).mean(axis=0)
    # related-key predictability
    worst = np.zeros(32)
    for bpos in range(80):
        d = x[:, bpos, :]
        for r in range(32):
            _, cnt = np.unique(d[:, r], return_counts=True)
            worst[r] = max(worst[r], cnt.max() / N11)
    # BLR
    a1, a2, a3 = (rand40(5000), rand40(5000)), (rand40(5000), rand40(5000)), (rand40(5000), rand40(5000))
    lhs = ksf(a1[0] ^ a2[0] ^ a3[0], a1[1] ^ a2[1] ^ a3[1]); rhs = ksf(*a1) ^ ksf(*a2) ^ ksf(*a3)
    # round-reduced cipher key avalanche
    _, r0 = present_encrypt(pt, base, record=True)
    _, r1 = present_encrypt(np.repeat(pt, 80), fl.reshape(-1, 32), record=True)
    chd = popcount(r1 ^ np.repeat(r0, 80, axis=0)) / 64
    E11[name] = {"ks_per_round_mean": hd.mean(axis=(0, 1)).tolist(), "ks_overall_mean": float(hd.mean()),
                 "sac_mean_abs_dev": float(np.abs(pm - 0.5).mean()),
                 "rk_max_prob_per_round": worst.tolist(), "rk_max_prob_rounds_2_32": float(worst[1:].max()),
                 "blr_full_linear_fraction": float(np.mean(np.all(lhs == rhs, axis=1))),
                 "blr_bits_satisfying": float((1 - popcount(lhs ^ rhs).sum(axis=1) / (64 * 32)).mean()),
                 "cipher_key_avalanche_per_round": chd.mean(axis=0).tolist()}
    print("E11", name, E11[name]["ks_overall_mean"], E11[name]["rk_max_prob_rounds_2_32"], np.round(chd.mean(axis=0)[:8], 3))
E11["PRESENT-80"]["sboxes"] = 31
E11["GAKS-PRESENT"]["sboxes"] = PRMP.sbox_count()
RES["E11"] = E11

# E12: the original paper's single-key bit-balance experiment, reproduced
E12 = {}
for name, ksf in KS.items():
    s = ksf(as_u64(0x3030303030303030))
    E12[name] = {"key_ascii_00000000_ones_per_subkey": popcount(s[0]).tolist()}
RES["E12"] = E12

json.dump(RES, open("results/eval.json", "w"), indent=1)
print("saved results/eval.json")

# E14: raw GA output (w=6) vs final design (w=7) on the evaluation keys
hd, p = ks_diffusion(lambda kk: ks_gaks(kk, PRM_GA), N2, FLIPS, 48, 16, des_keygen)
RES["E14_ga_raw_w6"] = {"w": PRM_GA.w, "sboxes": PRM_GA.sbox_count(),
    "significant_cells_per_round_bonferroni_1pct": (np.abs(p - 0.5) * 2 * np.sqrt(N2) > ZB).sum(axis=(0, 2)).tolist(),
    "sac_per_round_max": np.abs(p - 0.5).max(axis=(0, 2)).tolist()}
# E13: bit-level dependence between round keys of the same master key (max |corr| over bit pairs of different rounds)
def interkey_bitcorr(rk, nbits):
    from ga_design import bits_of as _b
    B = _b(rk, nbits).reshape(rk.shape[0], -1).astype(np.float32)
    B -= B.mean(0); B /= B.std(0) + 1e-9
    C = (B.T @ B) / B.shape[0]
    R_ = rk.shape[1]; blk = np.kron(np.eye(R_), np.ones((nbits, nbits))).astype(bool)
    C[blk] = 0
    return float(np.abs(C).max())
kk13 = rand_keys(40000)
RES["E13"] = {"N": 40000, "noise_ref": None,
              "DES": interkey_bitcorr(ks_des(kk13), 48), "GAKS-DES": interkey_bitcorr(KS["GAKS-DES"](kk13), 48)}
h13, l13 = rand40(20000), rand40(20000)
RES["E13"]["PRESENT-80 (K1..K8)"] = interkey_bitcorr(ks_present(h13, l13)[:, :8], 64)
RES["E13"]["GAKS-PRESENT (K1..K8)"] = interkey_bitcorr(ks_gaks_present(h13, l13, PRMP)[:, :8], 64)
RES["E13"]["noise_ref_DES_N40000"] = interkey_bitcorr(rng.integers(0, 2**48, (40000, 16), dtype=np.uint64), 48)
print("E13", RES["E13"]); print("E14", RES.get("E14_ga_raw_w6"))
json.dump(RES, open("results/eval.json", "w"), indent=1)
