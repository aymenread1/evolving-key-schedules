"""gift_sat.py -- SAT models for the minimum number of active S-boxes in r-round related-key differential
characteristics of GIFT-64, with the original key schedule and with GAKS-GIFT.

Each S-box (in the cipher and in the GAKS key schedule) is modelled exactly at the level of possible transitions:
one clause per impossible (input, output) difference pair of the GIFT difference distribution table, plus an
activity variable. Probabilities are not modelled, so every optimum is a valid lower bound on the number of
active S-boxes, and (6/16)^N bounds the probability of any characteristic with N active S-boxes.
The cardinality constraint uses an incremental totalizer; the bound k is raised until the instance is satisfiable.

Usage: python3 gift_sat.py <mode> <rounds> [time_limit_s]
  mode = sk-orig | rk-orig | rk-gaks | rk-gaks-free | ks-gaks
  rk-gaks counts every key-schedule S-box from the master key on (including the w warm-up updates);
  rk-gaks-free starts from an arbitrary non-zero key-state difference at the point where K1 is read.
Results are appended to results/gift/sat_rk.json."""
import json, os, sys, time, threading
import numpy as np
from pysat.solvers import Solver
from pysat.card import ITotalizer
from gift_core import GS, P64, GAKSGiftParams

S = [int(x) for x in GS]
DDT = np.zeros((16, 16), dtype=int)
for x in range(16):
    for a in range(16):
        DDT[a, S[x] ^ S[x ^ a]] += 1
IMPOSSIBLE = [(a, b) for a in range(16) for b in range(16) if DDT[a, b] == 0]
OUTF = os.environ.get("OUT", "results/gift/sat_rk.json")
OFFSET = int(os.environ.get("OFFSET", 0))   # rk-gaks-free: start the block at round 1+OFFSET
MAXP = DDT[1:, :].max() / 16                       # 6/16 for the GIFT S-box


class Model:
    def __init__(self):
        self.nv = 0; self.cl = []; self.act_cipher = []; self.act_ks = []; self.trace = []

    def var(self, n=1):
        v = list(range(self.nv + 1, self.nv + n + 1)); self.nv += n; return v

    def xor(self, a, b):                           # None = constant zero difference
        if a is None: return b
        if b is None: return a
        c, = self.var()
        self.cl += [[-a, -b, -c], [a, b, -c], [a, -b, c], [-a, b, c]]
        return c

    def sbox(self, xs, ks=False):
        if all(x is None for x in xs):
            return [None] * 4
        xs = [x if x is not None else self._zero() for x in xs]
        ys = self.var(4); A, = self.var()
        for v in xs + ys: self.cl.append([-v, A])
        self.cl.append([-A] + xs); self.cl.append([-A] + ys)
        for a, b in IMPOSSIBLE:
            pat = [(a >> q) & 1 for q in range(4)] + [(b >> q) & 1 for q in range(4)]
            self.cl.append([(-v if p else v) for v, p in zip(xs + ys, pat)])
        (self.act_ks if ks else self.act_cipher).append(A); self.trace.append(('sbox', ks, xs, ys))
        return ys

    def _zero(self):
        z, = self.var(); self.cl.append([-z]); return z


def gift_key_update(k):                            # k[i] = bit i of k7..k0 (bit 0 = LSB of k0)
    w = [k[16 * i:16 * i + 16] for i in range(8)]
    rot = lambda wd, s: [wd[(j + s) % 16] for j in range(16)]   # right rotation by s
    nw = w[2:8] + [rot(w[0], 12), rot(w[1], 2)]
    return [b for wd in nw for b in wd]


def gaks_step(m, st, t, prm, count=True):
    st = gift_key_update(st)
    H, L = st[64:], st[:64]; g = t % prm.nblocks
    for i in range(64):
        if (prm.masks[g] >> i) & 1: H[i], L[i] = L[i], H[i]
    L = [m.xor(L[i], H[(i - prm.alpha[g]) % 64]) for i in range(64)]
    H = [m.xor(H[i], L[(i - prm.beta[g]) % 64]) for i in range(64)]
    for half, sel in ((H, prm.smask >> 16), (L, prm.smask & 0xFFFF)):
        for j in range(16):
            if (sel >> j) & 1:
                half[4 * j:4 * j + 4] = m.sbox(half[4 * j:4 * j + 4], ks=True) if count else m.sbox(half[4 * j:4 * j + 4], ks=True)
    st = L + H; m.trace.append(('ks_state', t, st)); return st


def build(mode, r, prm=None):
    m = Model()
    if mode == "ks-gaks":                          # key schedule alone, r updates from the master key
        k = m.var(128); m.cl.append(list(k)); st = list(k)
        for t in range(r): st = gaks_step(m, st, t, prm)
        return m, k, None
    p = m.var(64); s = list(p)
    if mode == "sk-orig":
        k = None; m.cl.append(list(p))
    else:
        k = m.var(128); m.cl.append(list(k))         # related-key: non-zero key difference
    st = list(k) if k is not None else None
    t = 0
    if mode == "rk-gaks":
        for t in range(prm.w): st = gaks_step(m, st, t, prm)
        t = prm.w
    elif mode == "rk-gaks-free":
        t = prm.w + OFFSET                         # st is the key-state difference where K_{1+OFFSET} is read
    for rd in range(r):
        out = [None] * 64
        m.trace.append(('round_in', rd, list(s), list(st) if st is not None else None))
        for j in range(16):
            out[4 * j:4 * j + 4] = m.sbox(s[4 * j:4 * j + 4])
        m.trace.append(('round_out', rd, list(out)))
        nxt = [None] * 64
        for i in range(64): nxt[P64[i]] = out[i]
        if st is not None:
            for i in range(16):
                nxt[4 * i + 1] = m.xor(nxt[4 * i + 1], st[16 + i])   # U = k1
                nxt[4 * i] = m.xor(nxt[4 * i], st[i])                # V = k0
            if rd < r - 1:
                if mode == "rk-orig": st = gift_key_update(st)
                else: st = gaks_step(m, st, t, prm); t += 1
        s = nxt
    m.trace.append(('final', r, list(s)))
    return m, (p, k), s


def minimise(m, lb, tl, ub=80, progress=None):
    """Raise k from lb until SAT. CaDiCaL has no interrupt, so the time limit is enforced by the caller
    (each UNSAT answer is reported through `progress`, which records k+1 as a proven lower bound)."""
    acts = m.act_cipher + m.act_ks
    if not acts:
        return dict(optimal=True, value=0, seconds=0.0, model=None)
    tot = ITotalizer(lits=acts, ubound=min(ub, len(acts)), top_id=m.nv)
    sol = Solver(name="cadical195", bootstrap_with=m.cl + tot.cnf.clauses)
    t0 = time.time(); k = lb; res = None
    while k <= min(ub, len(acts)):
        ok = sol.solve(assumptions=[-tot.rhs[k]] if k < len(tot.rhs) else [])
        if ok:
            res = dict(optimal=True, value=k, seconds=round(time.time() - t0, 1), model=sol.get_model()); break
        k += 1
        if progress: progress(dict(optimal=False, value=None, lower_bound=k, seconds=round(time.time() - t0, 1)))
    if res is None:
        res = dict(optimal=False, value=None, lower_bound=k, seconds=round(time.time() - t0, 1), model=None)
    sol.delete()
    return res


def params():
    ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
    prm = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); prm.w = ev["final_w"]; return prm


if __name__ == "__main__":
    mode, rounds = sys.argv[1], [int(x) for x in sys.argv[2].split(",")]
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    lb0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    prm = params()
    try: out = json.load(open(OUTF))
    except Exception: out = {"solver": "CaDiCaL 1.9.5 via PySAT", "max_ddt_prob": MAXP, "rows": []}
    lb = lb0
    def save(row):
        out["rows"] = [q for q in out["rows"] if not (q["mode"] == row["mode"] and q["rounds"] == row["rounds"])] + [row]
        json.dump(out, open(OUTF, "w"), indent=1)
    for r in rounds:
        m, io, _ = build(mode, r, prm)
        base = dict(mode=mode + (f"@{OFFSET}" if mode == "rk-gaks-free" and OFFSET else ""), rounds=r, vars=m.nv, clauses=len(m.cl), sboxes_cipher=len(m.act_cipher), sboxes_ks=len(m.act_ks))
        save({**base, "optimal": False, "lower_bound": lb, "seconds": 0.0})
        res = minimise(m, lb, tl, progress=lambda p: (save({**base, **p}), print("  proven >=", p["lower_bound"], p["seconds"], "s", flush=True)))
        row = {**base, **{k: v for k, v in res.items() if k != "model"}}
        if res["model"] is not None:
            pos = set(l for l in res["model"] if l > 0)
            row["active_cipher"] = sum(a in pos for a in m.act_cipher); row["active_ks"] = sum(a in pos for a in m.act_ks)
            if mode != "ks-gaks":
                if io[1] is not None: row["delta_key"] = f"{sum(1 << i for i, v in enumerate(io[1]) if v in pos):032x}"
                row["delta_p"] = f"{sum(1 << i for i, v in enumerate(io[0]) if v in pos):016x}"
            else:
                row["delta_key"] = f"{sum(1 << i for i, v in enumerate(io) if v in pos):032x}"
        print(row, flush=True); save(row)
        lb = row["value"] if row.get("optimal") else row.get("lower_bound", lb)
