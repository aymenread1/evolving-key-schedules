"""
gift_design.py -- design-time GA for GAKS-GIFT (same fitness as ga_design.py, adapted to the 128-bit key state
and 28 x 32-bit round keys of GIFT-64). Usage: python3 gift_design.py  (env: NBLOCKS, LAMBDA, SUFFIX, GENS, POP)
"""
import json, os, time
import numpy as np
from gaks_core import U, popcount
from gift_core import *
from ga_design import bits_of

NB = int(os.environ.get("NBLOCKS", 2))
LAMBDA = float(os.environ.get("LAMBDA", 0.5)); CHI_TOL = 0.15; W_DET = 5.0; W_RHO = 5.0
POP, GENS = int(os.environ.get("POP", 40)), int(os.environ.get("GENS", 60))
RESUME = os.environ.get("RESUME")                 # path of a population checkpoint to continue from
TOUR, PC, ELITE, N = 3, 0.9, 2, 64
SEED = int(os.environ.get("SEED", 2026))
rk = np.random.default_rng(1)
DH, DL = rk.integers(0, 2**63, N, dtype=np.uint64) * U(2), rk.integers(0, 2**63, N, dtype=np.uint64) * U(2)
DH |= rk.integers(0, 2, N, dtype=np.uint64); DL |= rk.integers(0, 2, N, dtype=np.uint64)
NC = 16384
CH, CL = rk.integers(0, 2**63, NC, dtype=np.uint64) * U(2), rk.integers(0, 2**63, NC, dtype=np.uint64) * U(2)

def flips(H, L):
    n = len(H)
    FH = np.repeat(H[:, None], 128, 1).copy(); FL = np.repeat(L[:, None], 128, 1).copy()
    for b in range(64):
        FL[:, b] ^= U(1) << U(b); FH[:, 64 + b] ^= U(1) << U(b)
    return FH.ravel(), FL.ravel()
FH, FL = flips(DH, DL)

def interkey_corr(rk):
    w = popcount(rk).astype(np.float64); c = np.corrcoef(w.T) - np.eye(w.shape[1]); return float(np.nanmax(np.abs(c)))

def fitness(prm):
    base = ks_gaks_gift(DH, DL, prm)
    fl = ks_gaks_gift(FH, FL, prm).reshape(N, 128, ROUNDS)
    x = fl ^ base[:, None, :]
    a_r = (popcount(x) / 32).mean(axis=(0, 1))
    p = bits_of(x, 32).mean(axis=0); dev = np.abs(p - 0.5)
    sac = dev.mean() + dev.mean(axis=(0, 2)).max()
    det = (dev == 0.5).mean(axis=(0, 2)).max()
    chi = (4 * N * (p - 0.5) ** 2).mean(axis=(0, 2)) - 1.0
    ic = interkey_corr(ks_gaks_gift(CH, CL, prm))
    d = np.abs(a_r - 0.5)
    pen = d.mean() + d.max() + sac + W_DET * det + max(0.0, float(chi.max()) - CHI_TOL) + W_RHO * max(0.0, ic - 4 / np.sqrt(NC))
    cost = prm.sbox_count() / ((GAKSGiftParams.MAXW + ROUNDS - 1) * 32)
    return -pen - LAMBDA * cost, dict(a_min=float(a_r.min()), sac=float(sac), chi=float(chi.max()), interkey=ic, cost=prm.sbox_count())

def run():
    L = GAKSGiftParams.nbits(NB); PM = 2.0 / L; TAIL = 35; rng = np.random.default_rng(SEED)
    pop = rng.integers(0, 2, (POP, L), dtype=np.uint8); cache = {}
    g0 = 0
    if RESUME and os.path.exists(RESUME):
        ck = np.load(RESUME, allow_pickle=True).item(); pop = ck["pop"]; g0 = ck["gen"]; rng = np.random.default_rng(SEED + g0)
    F = lambda ch: cache.setdefault(ch.tobytes(), None) or cache.__setitem__(ch.tobytes(), fitness(GAKSGiftParams.from_bits(ch, NB))) or cache[ch.tobytes()]
    init = [F(c)[0] for c in pop]; hist = []; t0 = time.time()
    if RESUME and os.path.exists(RESUME):
        init = ck["init"]; hist = ck["hist"]
    for g in range(g0, GENS):
        sc = np.array([F(c)[0] for c in pop]); o = np.argsort(-sc)
        hist.append(dict(gen=g, best=float(sc[o[0]]), mean=float(sc.mean()), **F(pop[o[0]])[1]))
        print(f"gen {g} best {sc[o[0]]:.4f} {F(pop[o[0]])[1]} {time.time()-t0:.0f}s", flush=True)
        new = [pop[i].copy() for i in o[:ELITE]]
        while len(new) < POP:
            def tour():
                idx = rng.choice(POP, TOUR, replace=False); return pop[idx[np.argmax(sc[idx])]]
            a, b = tour().copy(), tour().copy()
            if rng.random() < PC:
                cp = rng.integers(1, L); a[cp:], b[cp:] = b[cp:].copy(), a[cp:].copy()
            for c in (a, b):
                m = rng.random(L) < PM; m[-TAIL:] |= rng.random(TAIL) < 0.05; c[m] ^= 1; new.append(c)
        pop = np.array(new[:POP])
        if RESUME:
            np.save(RESUME, {"pop": pop, "gen": g + 1, "init": init, "hist": hist}, allow_pickle=True)
    sc = np.array([F(c)[0] for c in pop]); best = pop[np.argmax(sc)]; prm = GAKSGiftParams.from_bits(best, NB)
    out = dict(nblocks=NB, lambda_=LAMBDA, seed_used=SEED, pop=POP, gens=GENS, seed=SEED, n_design_keys=N, best_fitness=float(sc.max()),
               best_info=F(best)[1], params=prm.to_dict(), chromosome="".join(map(str, best.tolist())), history=hist,
               random_designs=dict(mean=float(np.mean(init)), best=float(np.max(init))))
    os.makedirs("results/gift", exist_ok=True)
    json.dump(out, open(f"results/gift/ga_gift{os.environ.get('SUFFIX', '')}.json", "w"), indent=1)
    print("saved", prm.to_dict())

if __name__ == "__main__":
    run()
