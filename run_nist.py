"""
run_nist.py -- NIST SP 800-22 battery on round-key streams.
Usage: python3 run_nist.py <config>   config in: DES-counter GAKS-DES-counter DES-random GAKS-DES-random
                                                PRESENT-counter GAKS-PRESENT-counter
Each config: 100 sequences x 1,000,000 bits; round keys of consecutive master keys are concatenated.
"""
import json, sys, time
import numpy as np
from scipy.special import gammaincc
from gaks_core import *
import nist

import os
NSEQ, NBITS = int(os.environ.get("NSEQ", 100)), 1_000_000
PRM = final_des_params()
gp = json.load(open("results/ga_present.json")); PRMP = GAKSPresentParams.from_bits(gp["chromosome"])
POS = sorted(EFFECTIVE_BITS)  # non-parity bit positions, ascending

def counter_to_des_key(c):
    k = np.zeros_like(c)
    for j, p in enumerate(POS):
        k |= ((c >> U(j)) & U(1)) << U(p)
    return k

def to_bits(words, nbits):
    b = np.unpackbits(words.astype(">u8").view(np.uint8).reshape(words.shape + (8,)), axis=-1)[..., 64 - nbits:]
    return b.reshape(-1)

def gift_stream(cfg, s):
    import gift_core as GC
    per = -(-NBITS // 896); rng = np.random.default_rng(5000 + s)
    if cfg.endswith("counter"):
        L = np.arange(s * per, (s + 1) * per, dtype=np.uint64); H = np.zeros_like(L)
    else:
        H = rng.integers(0, 2**63, per, dtype=np.uint64) * U(2) + rng.integers(0, 2, per, dtype=np.uint64)
        L = rng.integers(0, 2**63, per, dtype=np.uint64) * U(2) + rng.integers(0, 2, per, dtype=np.uint64)
    if cfg.startswith("GAKSGIFT"):
        ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
        prm = GC.GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); prm.w = ev["final_w"]
        rk = GC.ks_gaks_gift(H, L, prm)
    else:
        rk = GC.ks_gift(H, L)
    return to_bits(rk, 32)[:NBITS]

def stream(cfg, s):
    if "GIFT" in cfg:
        return gift_stream(cfg, s)
    rng = np.random.default_rng(1000 + s)
    if "PRESENT" in cfg:
        per = -(-NBITS // 2048)
        c = np.arange(s * per, (s + 1) * per, dtype=np.uint64)
        kh, kl = np.zeros_like(c), c
        rk = ks_present(kh, kl) if cfg == "PRESENT-counter" else ks_gaks_present(kh, kl, PRMP)
        return to_bits(rk, 64)[:NBITS]
    per = -(-NBITS // 768)
    if cfg.endswith("counter"):
        k = counter_to_des_key(np.arange(s * per, (s + 1) * per, dtype=np.uint64))
    else:
        k = rng.integers(0, 2**63, per, dtype=np.uint64) * U(2) + rng.integers(0, 2, per, dtype=np.uint64)
    if cfg.startswith("GAKSHW"):
        g = json.load(open("results/hw_variant/ga_des_hw2.json"))
        hw = GAKSParams.from_bits(g["chromosome"], nblocks=g["nblocks"], rc=RC_CTR_DES)
        sk = ks_gaks(k, hw)
    else:
        sk = ks_des(k) if cfg.startswith("DES") else ks_gaks(k, PRM)
    return to_bits(sk, 48)[:NBITS]

def main(cfg):
    P = np.zeros((NSEQ, len(nist.TESTS)))
    t = time.time()
    for s in range(NSEQ):
        P[s] = nist.run_all(stream(cfg, s))
    alpha = 0.01
    pmin = (1 - alpha) - 3 * np.sqrt(alpha * (1 - alpha) / NSEQ)
    rows = []
    for j, name in enumerate(nist.TESTS):
        p = P[:, j]
        prop = float(np.mean(p >= alpha))
        hist = np.histogram(p, bins=10, range=(0, 1))[0]
        chi = np.sum((hist - NSEQ / 10) ** 2 / (NSEQ / 10))
        unif = float(gammaincc(9 / 2, chi / 2))
        rows.append(dict(test=name, proportion=prop, uniformity_p=unif,
                         passed=bool(prop >= pmin and unif >= 1e-4)))
    out = dict(config=cfg, nseq=NSEQ, nbits=NBITS, min_proportion=float(pmin), tests=rows,
               passed_count=sum(r["passed"] for r in rows), seconds=time.time() - t)
    json.dump(out, open(f"results/nist_{cfg}" + ("" if NSEQ == 100 else f"_n{NSEQ}") + ".json", "w"), indent=1)
    print(cfg, out["passed_count"], "/", len(rows), [(r["test"], r["proportion"], round(r["uniformity_p"], 4)) for r in rows])

if __name__ == "__main__":
    for cfg in sys.argv[1:]:
        main(cfg)
