"""evaluate_variant.py -- evaluation of the hardware-friendly GAKS-DES variants (cyclic gene blocks,
counter round constants). Step 1 applies the same validation rule as the main design (validation seed 777):
add one warm-up step if significant SAC cells remain after K1 or K1 has more than 5 biased cells.
Step 2 evaluates on the evaluation seed 20261122. Writes results/hw_variant/eval.json."""
import json, sys
import numpy as np
from gaks_core import *
from ga_design import bits_of, interkey_corr, CORR_KEYS
from linear_bypass import t_lin_des
from scipy.stats import norm

FLIPS = np.array([U(1) << U(b) for b in EFFECTIVE_BITS], dtype=np.uint64)
ZB = float(norm.isf(0.01 / 2 / (56 * 16 * 48)))

def load(name):
    g = json.load(open(f"results/hw_variant/ga_des_{name}.json"))
    return GAKSParams.from_bits(g["chromosome"], nblocks=g["nblocks"], rc=RC_CTR_DES if g["rcmode"] == "ctr" else None), g

def sac_sig(prm, rng, N=1000):
    k = rng.integers(0, 2**63, N, dtype=np.uint64) * U(2) + rng.integers(0, 2, N, dtype=np.uint64)
    b = ks_gaks(k, prm); f = ks_gaks((k[:, None] ^ FLIPS[None, :]).ravel(), prm).reshape(N, 56, 16)
    x = f ^ b[:, None, :]
    p = bits_of(x, 48).mean(0)
    hd = popcount(x) / 48
    return (np.abs(p - 0.5) * 2 * np.sqrt(N) > ZB).sum(axis=(0, 2)).tolist(), float(hd.mean()), float(np.abs(p - .5).mean()), k, b

out = {}
for name in sys.argv[1:] or ["hw4", "hw2"]:
    prm, g = load(name)
    sig_val, *_ = sac_sig(prm, np.random.default_rng(777))
    w0 = prm.w
    if sum(sig_val[1:]) > 0 or sig_val[0] > 5:
        prm.w = min(prm.w + 1, GAKSParams.MAXW)
    rng = np.random.default_rng(20261122)
    sig, av, sacdev, k, base = sac_sig(prm, rng)
    # E4 related-key predictability (single-bit differences)
    worst = 0.0
    for f in FLIPS:
        d = ks_gaks(k ^ f, prm) ^ base
        for r in range(16):
            _, c = np.unique(d[:, r], return_counts=True); worst = max(worst, c.max() / len(k))
    # E5 complementation, E6 weak keys
    kk = rng.integers(0, 2**63, 10000, dtype=np.uint64) * U(2); pp = rng.integers(0, 2**63, 10000, dtype=np.uint64)
    comp = float(np.mean(des_encrypt(~pp, ks_gaks(~kk, prm)) == ~des_encrypt(pp, ks_gaks(kk, prm))))
    weak = [0x0101010101010101, 0xFEFEFEFEFEFEFEFE, 0xE0E0E0E0F1F1F1F1, 0x1F1F1F1F0E0E0E0E]
    pt = rng.integers(0, 2**63, 1000, dtype=np.uint64)
    weak_inv = [float(np.mean(des_encrypt(des_encrypt(pt, np.repeat(ks_gaks(as_u64(w), prm), 1000, 0)), np.repeat(ks_gaks(as_u64(w), prm), 1000, 0)) == pt)) for w in weak]
    # E7 round-reduced key avalanche
    k7 = k[:300]; p7 = rng.integers(0, 2**63, 300, dtype=np.uint64)
    _, r0 = des_encrypt(p7, ks_gaks(k7, prm), record=True)
    _, r1 = des_encrypt(np.repeat(p7, 56), ks_gaks((k7[:, None] ^ FLIPS[None, :]).ravel(), prm), record=True)
    cav = (popcount(r1 ^ np.repeat(r0, 56, axis=0)) / 64).mean(axis=0)
    out[name] = dict(nblocks=g["nblocks"], rcmode=g["rcmode"], ga_w=w0, final_w=prm.w, sboxes=prm.sbox_count(),
                     smask=f"{prm.smask:04x}", validation_sig=sig_val, sig_per_round=sig, ks_avalanche=av, sac_mean_abs_dev=sacdev,
                     rk_max_prob=worst, complementation=comp, weak_key_involution=weak_inv,
                     rounds_to_049=int(next(i + 1 for i, v in enumerate(cav) if v >= 0.49)), cipher_av=cav[:4].tolist(),
                     T_lin=t_lin_des(prm), rho=interkey_corr(ks_gaks(CORR_KEYS, prm)))
    print(name, out[name])
json.dump(out, open("results/hw_variant/eval.json", "w"), indent=1)
