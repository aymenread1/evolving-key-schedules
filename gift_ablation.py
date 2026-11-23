"""gift_ablation.py -- cost-matched ablation for GAKS-GIFT: random parameter sets with the same structure
(two cyclic gene blocks, k S-boxes per step, given warm-up count) versus the GA design. Fresh keys (seed 4243).
Writes results/gift/ablation.json."""
import json
import numpy as np
from scipy.stats import norm
from gaks_core import U, popcount
from gift_core import *
from ga_design import bits_of
import gift_design as GD

N = 200
rng = np.random.default_rng(4243)
H = rng.integers(0, 2**63, N, dtype=np.uint64) * U(2); L = rng.integers(0, 2**63, N, dtype=np.uint64) * U(2)
FH, FL = GD.flips(H, L)
ZTHR = float(norm.isf(0.01 / 2 / (128 * ROUNDS * 32))) * 0.5 / np.sqrt(N)

def evaluate(prm):
    b = ks_gaks_gift(H, L, prm); x = ks_gaks_gift(FH, FL, prm).reshape(N, 128, ROUNDS) ^ b[:, None, :]
    a = (popcount(x) / 32).mean(axis=(0, 1)); d = np.abs(bits_of(x, 32).mean(0) - 0.5)
    sig = (d > ZTHR).sum(axis=(0, 2))
    return dict(a_min=float(a.min()), sig_K1=int(sig[0]), sig_later=int(sig[1:].sum()), fitness=float(GD.fitness(prm)[0]), cost=prm.sbox_count())

ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
GA = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"])
out = {}
for w in (GA.w, ev["final_w"]):
    ref = GAKSGiftParams(GA.masks, GA.alpha, GA.beta, GA.smask, w, GA.nblocks)
    out[f"w{w}"] = {"ga": evaluate(ref)}
    for k in (4, 5, 6):
        res = []
        for i in range(40):
            gg = np.random.default_rng(8000 + 100 * k + i)
            p = GAKSGiftParams.from_bits(gg.integers(0, 2, GAKSGiftParams.nbits(2)), 2)
            nib = gg.choice(32, k, replace=False); p.smask = sum(1 << int(j) for j in nib); p.w = w
            res.append(evaluate(p))
        ok = [r for r in res if r["a_min"] >= 0.49 and r["sig_later"] == 0]
        kga = bin(GA.smask).count("1")
        out[f"w{w}"][str(k)] = dict(cost=res[0]["cost"], frac_meets=len(ok) / len(res),
                                     median_sig_K1=float(np.median([r["sig_K1"] for r in res])),
                                     best_fitness=max(r["fitness"] for r in res), median_fitness=float(np.median([r["fitness"] for r in res])),
                                     n_better_than_ga=int(sum(r["fitness"] > out[f"w{w}"]["ga"]["fitness"] for r in res)) if k == kga else None)
        print(w, k, out[f"w{w}"][str(k)], flush=True)
    print(w, "GA", out[f"w{w}"]["ga"])
json.dump(out, open("results/gift/ablation.json", "w"), indent=1)
