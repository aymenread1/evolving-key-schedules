"""gift_milp.py -- MILP lower bounds (bit-oriented, S-box modelled by activity + differential branch number)
for GIFT-64: single-key and related-key characteristics with the original key schedule, and a key-schedule-only
model of GAKS-GIFT that checks T_lin - 1. Same modelling as milp_rk.py. Usage: python3 gift_milp.py <rounds> <tl>"""
import json, sys, time
import numpy as np
from milp_rk import Model, _ZERO, zero
from gift_core import GS, P64, GAKSGiftParams, ROUNDS

S = [int(x) for x in GS]
BN = min(bin(a).count("1") + bin(S[x] ^ S[x ^ a]).count("1") for a in range(1, 16) for x in range(16))

def impossible_transitions(box):
    ddt = np.zeros((16, 16), dtype=int)
    for x in range(16):
        for a in range(16): ddt[a, box[x] ^ box[x ^ a]] += 1
    return [(a, b) for a in range(16) for b in range(16) if ddt[a, b] == 0]
IMPOSSIBLE = impossible_transitions(S)

class M(Model):
    """S-box modelled exactly at the level of possible transitions: every impossible (dx, dy) pair of the
    difference distribution table is excluded by one no-good cut (probabilities are not modelled)."""
    def sbox(self, xs, cost=True):
        ys = self.var(4); A, = self.var(1)
        for v in xs + ys: self.con({A: 1, v: -1}, 0, np.inf)
        self.con({**{v: 1 for v in xs}, A: -1}, 0, np.inf); self.con({**{v: 1 for v in ys}, A: -1}, 0, np.inf)
        self.con({**{v: 1 for v in xs + ys}, A: -BN}, 0, np.inf)
        bitsv = xs + ys                                   # xs[q], ys[q] = bit q (LSB first)
        for a, b in IMPOSSIBLE:
            pat = [(a >> q) & 1 for q in range(4)] + [(b >> q) & 1 for q in range(4)]
            coef, ones = {}, 0
            for v, pb in zip(bitsv, pat):
                if pb: coef[v] = coef.get(v, 0) - 1; ones += 1
                else: coef[v] = coef.get(v, 0) + 1
            self.con(coef, 1 - ones, np.inf)              # pattern != (a, b)
        if cost: self.obj[A] = self.obj.get(A, 0) + 1
        return ys, A

def gift_key_update(k):     # k: list of 128 entries, index = bit position of k7..k0 (bit 0 = LSB of k0)
    w = [k[16 * i:16 * i + 16] for i in range(8)]
    rot = lambda wd, s: [wd[(j + s) % 16] for j in range(16)]   # right rotation by s: new bit j = old bit j+s
    nw = w[2:8] + [rot(w[0], 12), rot(w[1], 2)]
    return [b for wd in nw for b in wd]

def build(r, related):
    _ZERO.clear(); m = M()
    p = m.var(64); k = m.var(128) if related else [None] * 128
    s = list(p); act = []
    for rd in range(r):
        out = [None] * 64
        for j in range(16):
            xs = [s[4 * j + q] if s[4 * j + q] is not None else zero(m) for q in range(4)]
            ys, A = m.sbox(xs); act.append(A); out[4 * j:4 * j + 4] = ys
        t = [None] * 64
        for i in range(64): t[P64[i]] = out[i]
        if related:
            for i in range(16):
                for pos, kb in ((4 * i + 1, k[16 + i]), (4 * i, k[i])):   # U = k1 -> b(4i+1), V = k0 -> b(4i)
                    t[pos] = m.xor(t[pos], kb)
            k = gift_key_update(k)
        s = t
    m.con({v: 1 for v in p + ([v for v in k if v is not None] if related else [])}, 1, np.inf)
    if related:
        pass
    return m, act

def solve(m, tl):
    t0 = time.time(); res = m.solve(tl)
    return dict(optimal=bool(res.status == 0), value=float(res.fun) if res.x is not None else None, seconds=round(time.time() - t0, 1))

def ks_only_gaks(prm, steps, tl):
    _ZERO.clear(); m = M(); k = m.var(128); m.con({v: 1 for v in k}, 1, np.inf); st = list(k)
    for t in range(steps):
        st = gift_key_update(st)
        H, L = st[64:], st[:64]; g = t % prm.nblocks
        for i in range(64):
            if (prm.masks[g] >> i) & 1: H[i], L[i] = L[i], H[i]
        L = [m.xor(L[i], H[(i - prm.alpha[g]) % 64]) for i in range(64)]
        H = [m.xor(H[i], L[(i - prm.beta[g]) % 64]) for i in range(64)]
        for half, sel in ((H, prm.smask >> 16), (L, prm.smask & 0xFFFF)):
            for j in range(16):
                if (sel >> j) & 1:
                    xs = [x if x is not None else zero(m) for x in half[4 * j:4 * j + 4]]
                    ys, A = m.sbox(xs); half[4 * j:4 * j + 4] = ys
        st = L + H
    return solve(m, tl)

if __name__ == "__main__":
    rounds = [int(x) for x in sys.argv[1].split(",")]; tl = float(sys.argv[2])
    try: out = json.load(open("results/gift/milp.json"))
    except Exception: out = {"branch_number": BN, "rows": []}
    for r in rounds:
        for rel in (False, True):
            m, act = build(r, rel); row = dict(rounds=r, related_key=rel, **solve(m, tl))
            out["rows"] = [q for q in out["rows"] if not (q["rounds"] == r and q["related_key"] == rel)] + [row]
            print(row, flush=True); json.dump(out, open("results/gift/milp.json", "w"), indent=1)
    if len(sys.argv) > 3:
        ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
        prm = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); prm.w = ev["final_w"]
        for st in [int(x) for x in sys.argv[3].split(",")]:
            row = dict(gaks_ks_only_steps=st, **ks_only_gaks(prm, st, tl)); out.setdefault("gaks_ks_only", []).append(row)
            print(row, flush=True); json.dump(out, open("results/gift/milp.json", "w"), indent=1)
