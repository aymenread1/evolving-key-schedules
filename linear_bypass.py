"""
linear_bypass.py -- exact GF(2) analysis of S-box-free key-difference trails in key schedules.

A key difference that never reaches an S-box input propagates through the key schedule linearly and
with probability 1. The set of such differences surviving t state updates is the kernel of a linear
system: "all S-box input bits at steps 1..t have zero difference". T_lin is the first t at which this
kernel is {0}, i.e. after which every non-zero key difference has activated at least one key-schedule
S-box. Lower bound: T_lin >= ceil(n / (4k)) for an n-bit state and k S-boxes per step.
Writes results/linear_bypass.json.
"""
import json, math
import numpy as np
from gaks_core import *

def rank_gf2(rows):
    rows = [r for r in rows if r]; rank = 0; piv = {}
    for r in rows:
        while r:
            h = r.bit_length() - 1
            if h in piv: r ^= piv[h]
            else: piv[h] = r; rank += 1; break
    return rank

def basis(n):
    return [1 << i for i in range(n)]

def t_lin_des(prm, max_steps=40, start=0):
    """State = 56-bit (C<<28)|D after PC-1; basis vector i = difference in state bit i."""
    n = 56
    C = np.array([(1 << i) >> 28 for i in range(n)], dtype=object)
    D = np.array([(1 << i) & ((1 << 28) - 1) for i in range(n)], dtype=object)
    rot = lambda x, s: ((x << s) | (x >> (28 - s))) & ((1 << 28) - 1)
    rows = []
    nib = [j for j in range(14) if (prm.smask >> j) & 1]
    for t in range(start, min(start + max_steps, prm.w + 16)):
        r = t - prm.w; s = SHIFTS[r] if r >= 0 else 1
        C = np.array([rot(int(c), s) for c in C], dtype=object); D = np.array([rot(int(d), s) for d in D], dtype=object)
        m = prm.masks[t]; d = (C ^ D) & m; C, D = C ^ d, D ^ d
        D = np.array([int(dd) ^ rot(int(cc), prm.alpha[t]) for cc, dd in zip(C, D)], dtype=object)
        C = np.array([int(cc) ^ rot(int(dd), prm.beta[t]) for cc, dd in zip(C, D)], dtype=object)
        S = [(int(c) << 28) | int(dd) for c, dd in zip(C, D)]
        for j in nib:
            for q in range(4):
                b = 4 * j + q
                rows.append(sum(((S[i] >> b) & 1) << i for i in range(n)))
        if rank_gf2(rows) == n:
            return t + 1 - start
    return None

def t_lin_present(kind, prm=None, max_steps=60, start=0, total=None):
    n = 80; M = (1 << 40) - 1
    H = [((1 << i) >> 40) & M for i in range(n)]; L = [(1 << i) & M for i in range(n)]
    def rot80(h, l, k):
        v = (h << 40) | l; v = ((v << k) | (v >> (80 - k))) & ((1 << 80) - 1); return v >> 40, v & M
    rot40 = lambda x, s: ((x << s) | (x >> (40 - s))) & M
    rows = []
    for t in range(start, start + max_steps if total is None else total):
        HL = [rot80(h, l, 61) for h, l in zip(H, L)]; H = [a for a, _ in HL]; L = [b for _, b in HL]
        if kind == "present":
            bits = [76, 77, 78, 79]
        else:
            g = t % prm.NMASK; m = prm.masks[g]
            nh, nl = [], []
            for h, l in zip(H, L):
                d = (h ^ l) & m; h ^= d; l ^= d
                l ^= rot40(h, prm.alpha[g]); h ^= rot40(l, prm.beta[g]); nh.append(h); nl.append(l)
            H, L = nh, nl
            bits = [40 + 4 * j + q for j in range(10) if (prm.smask >> (10 + j)) & 1 for q in range(4)] + \
                   [4 * j + q for j in range(10) if (prm.smask >> j) & 1 for q in range(4)]
        V = [(h << 40) | l for h, l in zip(H, L)]
        for b in bits:
            rows.append(sum(((V[i] >> b) & 1) << i for i in range(n)))
        if rank_gf2(rows) == n:
            return t + 1 - start
    return None

def guaranteed_active(tfun, total):
    """Greedy tiling: windows [s, s+T(s)) each contain >= 1 active key-schedule S-box."""
    s, cnt, windows = 0, 0, []
    while True:
        T = tfun(s)
        if T is None or s + T > total:
            return cnt, windows
        cnt += 1; windows.append(T); s += T

if __name__ == "__main__":
    P = final_des_params()
    GP = GAKSPresentParams.from_bits(json.load(open("results/ga_present.json"))["chromosome"])
    kd = bin(P.smask).count("1"); kp = bin(GP.smask).count("1")
    out = {"GAKS-DES": {"sboxes_per_step": kd, "T_lin": t_lin_des(P), "lower_bound": math.ceil(56 / (4 * kd)),
                        "steps_before_K1": P.w + 1, "steps_total": P.w + 16},
           "GAKS-PRESENT": {"sboxes_per_step": kp, "T_lin": t_lin_present("gaks", GP), "lower_bound": math.ceil(80 / (4 * kp)),
                            "steps_before_K1": GP.w, "steps_total": GP.w + 31},
           "PRESENT-80": {"sboxes_per_step": 1, "T_lin": t_lin_present("present"), "lower_bound": 20, "steps_total": 31},
           "DES": {"sboxes_per_step": 0, "T_lin": None, "note": "no S-box: every key difference propagates linearly"}}
    gd, wd = guaranteed_active(lambda st: t_lin_des(P, start=st), P.w + 16)
    gp, wp = guaranteed_active(lambda st: t_lin_present("gaks", GP, start=st, total=GP.w + 31), GP.w + 31)
    pp, wpp = guaranteed_active(lambda st: t_lin_present("present", start=st, total=31), 31)
    out["GAKS-DES"].update(T_lin_all_starts=[t_lin_des(P, start=st) for st in range(P.w + 16)], guaranteed_active=gd, windows=wd)
    out["GAKS-PRESENT"].update(guaranteed_active=gp, windows=wp)
    out["PRESENT-80"].update(guaranteed_active=pp, windows=wpp)
    # trade-off: random GAKS-DES designs with k S-boxes per step (w = 7)
    trade = {}
    for k in (1, 2, 3, 4, 5, 7, 14):
        vals = []
        for i in range(10):
            g = np.random.default_rng(7000 + 100 * k + i)
            bits = g.integers(0, 2, GAKSParams.NBITS)
            q = GAKSParams.from_bits(bits); nib = g.choice(14, k, replace=False)
            q.smask = sum(1 << int(j) for j in nib); q.w = 7
            vals.append(t_lin_des(q))
        trade[k] = {"cost": 23 * k, "bound": math.ceil(56 / (4 * k)), "T_lin_min": min(v for v in vals if v), "T_lin_median": float(np.median([v or 99 for v in vals]))}
    out["des_tradeoff_random_designs"] = trade
    json.dump(out, open("results/linear_bypass.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
