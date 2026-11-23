"""
ga_design.py -- Offline (design-time) genetic algorithm that selects the parameters of the
GAKS key schedule. The GA is run ONCE; its output is a fixed public parameter set, so the
resulting key schedule is deterministic. Usage: python3 ga_design.py [des|present]
"""
import json, sys, time
import numpy as np
from gaks_core import *

import os
LAMBDA = float(os.environ.get("LAMBDA", 0.50))      # weight of the implementation-cost penalty
CHI_TOL = float(os.environ.get("CHI_TOL", 0.15))    # tolerance of the per-round chi-square SAC term
W_DET = float(os.environ.get("W_DET", 5.0))         # weight of deterministic-cell term
W_RHO = float(os.environ.get("W_RHO", 5.0))         # weight of round-key correlation term
POP, GENS = 50, 80
TOUR, PC, ELITE = 3, 0.9, 2
N_DESIGN_KEYS = 128
SEED = 2026

rng_keys = np.random.default_rng(1)          # design key set (disjoint from evaluation seeds)
DES_KEYS = rng_keys.integers(0, 2**63, N_DESIGN_KEYS, dtype=np.uint64) * U(2) + rng_keys.integers(0, 2, N_DESIGN_KEYS, dtype=np.uint64)
FLIP_DES = [U(1) << U(b) for b in EFFECTIVE_BITS]
NC = 16384                                     # keys for the inter-round-key correlation term
CORR_KEYS = rng_keys.integers(0, 2**63, NC, dtype=np.uint64) * U(2)
CORR_H = rng_keys.integers(0, 2**40, NC, dtype=np.uint64); CORR_L = rng_keys.integers(0, 2**40, NC, dtype=np.uint64)

def interkey_corr(rk):
    """max |corr| between Hamming weights of distinct round keys (dependence between round keys)."""
    w = popcount(rk).astype(np.float64)
    c = np.corrcoef(w.T) - np.eye(w.shape[1])
    return float(np.nanmax(np.abs(c)))
P_H = rng_keys.integers(0, 2**40, N_DESIGN_KEYS, dtype=np.uint64)
P_L = rng_keys.integers(0, 2**40, N_DESIGN_KEYS, dtype=np.uint64)

def bits_of(arr, nbits):
    """(..., ) uint64 -> (..., nbits) uint8, LSB-first is irrelevant for SAC statistics."""
    b = np.unpackbits(arr.astype(">u8").view(np.uint8).reshape(arr.shape + (8,)), axis=-1)
    return b[..., 64 - nbits:]

def diffusion_metrics(base, flipped, nbits):
    """base (N,R), flipped (N,F,R) -> per-round avalanche, SAC deviation."""
    x = flipped ^ base[:, None, :]
    aval = popcount(x) / nbits                           # (N,F,R)
    a_r = aval.mean(axis=(0, 1))                         # (R,)
    p = bits_of(x, nbits).mean(axis=0)                   # (F,R,nbits)
    dev = np.abs(p - 0.5)
    sac = dev.mean() + dev.mean(axis=(0, 2)).max()          # overall + worst round
    det = (dev == 0.5).mean(axis=(0, 2)).max()              # deterministic (linear) cells, worst round
    # per-round chi-square SAC statistic: mean of z^2 = 4N(p-1/2)^2 is ~1 for an ideal function
    chi = (4 * x.shape[0] * (p - 0.5) ** 2).mean(axis=(0, 2)) - 1.0
    chi_pen = max(0.0, float(chi.max()) - CHI_TOL)             # tolerance ~5 standard deviations
    return a_r, sac + W_DET * det + chi_pen

def des_fitness(prm, keys=DES_KEYS):
    base = ks_gaks(keys, prm)
    fk = (keys[:, None] ^ np.array(FLIP_DES, dtype=np.uint64)[None, :]).ravel()
    fl = ks_gaks(fk, prm).reshape(len(keys), 56, 16)
    a_r, sac = diffusion_metrics(base, fl, 48)
    dev = np.abs(a_r - 0.5)
    cost = prm.sbox_count() / ((GAKSParams.MAXW + 16) * 14)
    ic = interkey_corr(ks_gaks(CORR_KEYS, prm))
    icp = W_RHO * max(0.0, ic - 4.0 / np.sqrt(NC))
    return -(dev.mean() + dev.max() + sac + icp) - LAMBDA * cost, dict(a_min=float(a_r.min()), sac=float(sac), interkey=ic, cost=prm.sbox_count())

def present_fitness(prm, kh=P_H, kl=P_L):
    base = ks_gaks_present(kh, kl, prm)
    fh, flo = [], []
    for b in range(80):
        if b >= 40:
            fh.append(kh ^ U(1 << (b - 40))); flo.append(kl)
        else:
            fh.append(kh); flo.append(kl ^ U(1 << b))
    fh = np.stack(fh, 1).ravel(); flo = np.stack(flo, 1).ravel()
    fl = ks_gaks_present(fh, flo, prm).reshape(len(kh), 80, 32)
    a_r, sac = diffusion_metrics(base, fl, 64)
    dev = np.abs(a_r - 0.5)
    cost = prm.sbox_count() / ((GAKSPresentParams.MAXW + 31) * 20)
    ic = interkey_corr(ks_gaks_present(CORR_H, CORR_L, prm))
    icp = W_RHO * max(0.0, ic - 4.0 / np.sqrt(NC))
    return -(dev.mean() + dev.max() + sac + icp) - LAMBDA * cost, dict(a_min=float(a_r.min()), sac=float(sac), interkey=ic, cost=prm.sbox_count())

def run(kind):
    P = GAKSParams if kind == "des" else GAKSPresentParams
    NB = int(os.environ.get("NBLOCKS", "0")) or None             # cyclic gene blocks (DES only)
    RCMODE = os.environ.get("RCMODE", "sha")                     # "sha" or "ctr" (DES only)
    if kind == "des" and (NB or RCMODE != "sha"):
        class P:                                                 # wrapper with the restricted encoding
            NBITS = GAKSParams.nbits(NB)
            @staticmethod
            def from_bits(ch):
                return GAKSParams.from_bits(ch, nblocks=NB, rc=(RC_CTR_DES if RCMODE == "ctr" else None))
    fit = des_fitness if kind == "des" else present_fitness
    L = P.NBITS
    PM = 2.0 / L
    TAIL = 17 if kind == "des" else 23          # structural genes: S-box mask + warm-up count
    PM_TAIL = 0.05                              # higher mutation rate on structural genes
    rng = np.random.default_rng(SEED)
    pop = rng.integers(0, 2, (POP, L), dtype=np.uint8)
    cache = {}
    def F(ch):
        k = ch.tobytes()
        if k not in cache:
            cache[k] = fit(P.from_bits(ch))
        return cache[k]
    hist = []
    init_scores = [F(c)[0] for c in pop]          # = fitness of random (non-optimised) designs
    t0 = time.time()
    for g in range(GENS):
        scores = np.array([F(c)[0] for c in pop])
        order = np.argsort(-scores)
        best = pop[order[0]]
        hist.append(dict(gen=g, best=float(scores[order[0]]), mean=float(scores.mean()), **F(best)[1]))
        print(f"{kind} gen {g:2d} best {scores[order[0]]:.4f} mean {scores.mean():.4f} {F(best)[1]} {time.time()-t0:.0f}s", flush=True)
        new = [pop[i].copy() for i in order[:ELITE]]                    # elitism
        while len(new) < POP:
            def tour():                                                  # selection
                idx = rng.choice(POP, TOUR, replace=False)
                return pop[idx[np.argmax(scores[idx])]]
            a, b = tour().copy(), tour().copy()
            if rng.random() < PC:                                        # single-point crossover
                cp = rng.integers(1, L)
                a[cp:], b[cp:] = b[cp:].copy(), a[cp:].copy()
            for c in (a, b):                                             # bit-flip mutation
                m = rng.random(L) < PM
                m[-TAIL:] |= rng.random(TAIL) < PM_TAIL
                c[m] ^= 1
                new.append(c)
        pop = np.array(new[:POP])
    scores = np.array([F(c)[0] for c in pop])
    best = pop[np.argmax(scores)]
    prm = P.from_bits(best)
    out = dict(nblocks=NB if kind == "des" else None, rcmode=RCMODE if kind == "des" else None, kind=kind, lambda_=LAMBDA, chi_tol=CHI_TOL, w_det=W_DET, w_rho=W_RHO, pm_tail=PM_TAIL, pop=POP, gens=GENS, seed=SEED, n_design_keys=N_DESIGN_KEYS,
               best_fitness=float(scores.max()), best_info=F(best)[1], params=prm.to_dict(),
               chromosome="".join(map(str, best.tolist())), history=hist,
               random_designs=dict(n=len(init_scores), mean=float(np.mean(init_scores)), std=float(np.std(init_scores)),
                                   best=float(np.max(init_scores)), worst=float(np.min(init_scores))))
    suffix = os.environ.get("SUFFIX", "" if LAMBDA == 0.50 else f"_lambda{LAMBDA}")
    json.dump(out, open(f"results/ga_{kind}{suffix}.json", "w"), indent=1)
    print("saved", prm.to_dict())

if __name__ == "__main__":
    import os; os.makedirs("results", exist_ok=True)
    run(sys.argv[1] if len(sys.argv) > 1 else "des")
