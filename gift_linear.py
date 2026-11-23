"""gift_linear.py -- exact S-box bypass length T_lin (GF(2)) for GAKS-GIFT, as in linear_bypass.py.
Key state = 128-bit integer (H<<64)|L. Writes results/gift/linear_bypass.json."""
import json, math
import numpy as np
from gift_core import *
from linear_bypass import rank_gf2

M64i = (1 << 64) - 1; M16i = 0xFFFF
r16 = lambda x, s: ((x >> s) | (x << (16 - s))) & M16i
rl64 = lambda x, s: ((x << s) | (x >> (64 - s))) & M64i

def gift_update_int(v):
    H, L = v >> 64, v & M64i
    k0, k1 = L & M16i, (L >> 16) & M16i
    nL = (L >> 32) | ((H & 0xFFFFFFFF) << 32)
    nH = (H >> 32) | (r16(k0, 12) << 32) | (r16(k1, 2) << 48)
    return (nH << 64) | nL

def lin_step(v, t, prm):
    v = gift_update_int(v); H, L = v >> 64, v & M64i
    g = t % prm.nblocks; d = (H ^ L) & prm.masks[g]; H ^= d; L ^= d
    L ^= rl64(H, prm.alpha[g]); H ^= rl64(L, prm.beta[g])
    return (H << 64) | L

def sbox_bits(prm):
    hi = [64 + 4 * j + q for j in range(16) if (prm.smask >> (16 + j)) & 1 for q in range(4)]
    lo = [4 * j + q for j in range(16) if (prm.smask >> j) & 1 for q in range(4)]
    return hi + lo

def t_lin(prm, start, total):
    basis = [1 << i for i in range(128)]
    # bring the basis to the state before step `start` is irrelevant: the kernel is over differences entering step `start`
    rows, bits = [], sbox_bits(prm)
    for t in range(start, total):
        basis = [lin_step(v, t, prm) for v in basis]
        for b in bits:
            rows.append(sum(((basis[i] >> b) & 1) << i for i in range(128)))
        if rank_gf2(rows) == 128:
            return t + 1 - start
    return None

if __name__ == "__main__":
    ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
    prm = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); prm.w = ev["final_w"]
    total = prm.w + ROUNDS - 1
    k = bin(prm.smask).count("1")
    Ts = [t_lin(prm, s, total) for s in range(total)]
    cnt, s, win = 0, 0, []
    while True:
        T = Ts[s] if s < total else None
        if T is None or s + T > total: break
        cnt += 1; win.append(T); s += T
    out = {"sboxes_per_step": k, "bound": math.ceil(128 / (4 * k)), "T_lin": Ts[0], "T_lin_all_starts": Ts,
           "steps_total": total, "guaranteed_active": cnt, "windows": win,
           "GIFT-64_original": "no S-box in the key schedule: T_lin infinite"}
    json.dump(out, open("results/gift/linear_bypass.json", "w"), indent=1)
    print(out)
