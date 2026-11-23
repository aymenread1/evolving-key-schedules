"""
milp_rk.py -- MILP lower bounds on the number of active S-boxes in related-key differential
trails of r-round PRESENT-80, with the original key schedule and with GAKS-PRESENT.

Bit-oriented truncated model (Mouha et al., Inscrypt 2011): one binary variable per bit
difference; XOR modelled by the four facets of its convex hull; every S-box (cipher and key schedule) gets an
activity variable A with input-nonzero <=> output-nonzero and the differential branch number
of the PRESENT S-box (=3, computed below). This relaxes the real S-box, so the optimum is a
valid LOWER bound on the number of active S-boxes of any related-key characteristic.
Objective: cipher + key-schedule active S-boxes. Usage: python3 milp_rk.py <rmax> <time_limit_s>
"""
import json, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
from gaks_core import PRESENT_S, PLAYER, GAKSPresentParams

S = [int(x) for x in PRESENT_S]
BN = min(bin(a).count("1") + bin(S[x] ^ S[x ^ a]).count("1") for a in range(1, 16) for x in range(16))

class Model:
    def __init__(self):
        self.nv = 0; self.rows = []; self.lo = []; self.hi = []; self.obj = {}
    def var(self, n=1):
        v = list(range(self.nv, self.nv + n)); self.nv += n; return v
    def con(self, coef, lo, hi):
        self.rows.append(coef); self.lo.append(lo); self.hi.append(hi)
    def xor(self, a, b):
        """c = a xor b on difference bits, modelled by the four facets of its convex hull."""
        if a is None: return b
        if b is None: return a
        c, = self.var(1)                     # convex hull of {(a,b,c): c = a xor b}
        self.con({a: 1, b: 1, c: -1}, 0, np.inf)
        self.con({a: 1, b: -1, c: 1}, 0, np.inf)
        self.con({a: -1, b: 1, c: 1}, 0, np.inf)
        self.con({a: 1, b: 1, c: 1}, -np.inf, 2)
        return c
    def sbox(self, xs, cost=True):
        ys = self.var(4); A, = self.var(1)
        for v in xs + ys:
            self.con({A: 1, v: -1}, 0, np.inf)                  # A >= bit
        self.con({**{v: 1 for v in xs}, A: -1}, 0, np.inf)     # sum(in) >= A
        self.con({**{v: 1 for v in ys}, A: -1}, 0, np.inf)     # sum(out) >= A
        self.con({**{v: 1 for v in xs + ys}, A: -BN}, 0, np.inf)  # branch number
        self.obj[A] = self.obj.get(A, 0) + 1
        return ys, A
    def solve(self, tl):
        r, c, v = [], [], []
        for i, row in enumerate(self.rows):
            for j, a in row.items():
                r.append(i); c.append(j); v.append(a)
        Am = coo_matrix((v, (r, c)), shape=(len(self.rows), self.nv)).tocsr()
        cvec = np.zeros(self.nv)
        for j, a in self.obj.items(): cvec[j] = a
        res = milp(cvec, constraints=LinearConstraint(Am, self.lo, self.hi), integrality=np.ones(self.nv),
                   bounds=Bounds(0, 1), options={"time_limit": tl, "disp": False, "mip_rel_gap": 0})
        return res

# ---- key-schedule state as 80 entries (index = bit position of the 80-bit register, LSB = 0); None = zero difference
def rot80(st, n):
    return [st[(j - n) % 80] for j in range(80)]

def present_ks_step(m, st, i, sbox_active):
    st = rot80(st, 61)
    xs = [st[76 + q] for q in range(4)]
    if all(x is None for x in xs):
        return st
    xs = [x if x is not None else zero(m) for x in xs]
    ys, A = m.sbox(xs); sbox_active.append(A)
    st = st[:76] + ys
    return st

_ZERO = {}
def zero(m):
    if id(m) not in _ZERO:
        z, = m.var(1); m.con({z: 1}, 0, 0); _ZERO[id(m)] = z
    return _ZERO[id(m)]

def gaks_step(m, st, t, prm, sbox_active):
    st = rot80(st, 61)
    H, L = st[40:], st[:40]                       # H[i] = register bit 40+i
    g = t % prm.NMASK; mask = prm.masks[g]
    for i in range(40):
        if (mask >> i) & 1:
            H[i], L[i] = L[i], H[i]
    a, b = prm.alpha[g], prm.beta[g]
    L = [m.xor(L[i], H[(i - a) % 40]) for i in range(40)]
    H = [m.xor(H[i], L[(i - b) % 40]) for i in range(40)]
    for half, sel in ((H, prm.smask >> 10), (L, prm.smask & 0x3FF)):
        for j in range(10):
            if (sel >> j) & 1:
                xs = half[4 * j:4 * j + 4]
                if all(x is None for x in xs):
                    continue
                xs = [x if x is not None else zero(m) for x in xs]
                ys, A = m.sbox(xs); sbox_active.append(A)
                half[4 * j:4 * j + 4] = ys
    return L + H

def round_key(st):
    return st[16:80]                               # bits 79..16 -> round-key bits 63..0

def build(r, schedule, prm=None, single_key=False, key_diff=False):
    _ZERO.clear()
    m = Model()
    k0 = [None] * 80 if single_key else m.var(80)
    p0 = m.var(64)
    ks_act, cs_act = [], []
    # round keys K1..Kr
    st = list(k0); keys = []
    if schedule == "present":
        keys.append(round_key(st))
        for i in range(1, r):
            st = present_ks_step(m, st, i, ks_act); keys.append(round_key(st))
    else:
        if prm.w == 0: keys.append(round_key(st))
        t = 0
        while len(keys) < r:
            st = gaks_step(m, st, t, prm, ks_act)
            if t + 1 >= prm.w: keys.append(round_key(st))
            t += 1
    s = list(p0)
    for i in range(r):
        u = [m.xor(s[b], keys[i][b]) if keys[i][b] is not None else s[b] for b in range(64)]
        out = [None] * 64
        for j in range(16):
            xs = [u[4 * j + q] if u[4 * j + q] is not None else zero(m) for q in range(4)]
            ys, A = m.sbox(xs); cs_act.append(A)
            out[4 * j:4 * j + 4] = ys
        s = [None] * 64
        for b in range(64):
            s[PLAYER[b]] = out[b]
    nz = {v: 1 for v in p0 + ([] if single_key else [x for x in k0])}
    m.con(nz, 1, np.inf)                           # non-trivial difference
    if key_diff:
        m.con({v: 1 for v in k0}, 1, np.inf)        # force a non-zero master-key difference
    return m, ks_act, cs_act

CFGS = {"PRESENT-80 single-key": ("present", True, False),
        "PRESENT-80 related-key": ("present", False, True),
        "GAKS-PRESENT related-key": ("gaks", False, True)}

def run(rounds, tl, names):
    g = json.load(open("results/ga_present.json")); prm = GAKSPresentParams.from_bits(g["chromosome"])
    try:
        out = json.load(open("results/milp_rk.json"))
    except Exception:
        out = {"branch_number": BN, "rows": []}
    for r in rounds:
        for name in names:
            sch, sk, kd = CFGS[name]
            t0 = time.time()
            m, ks, cs = build(r, sch, prm, sk, kd)
            res = m.solve(tl)
            x = res.x
            row = dict(config=name, rounds=r, vars=m.nv, time_limit_s=tl, optimal=bool(res.status == 0),
                       seconds=round(time.time() - t0, 1),
                       best_found=float(res.fun) if x is not None else None,
                       lower_bound=float(res.fun) if res.status == 0 else float(getattr(res, "mip_dual_bound", np.nan)))
            if x is not None:
                row["cipher_active"] = int(round(sum(x[a] for a in cs))); row["ks_active"] = int(round(sum(x[a] for a in ks)))
            out["rows"] = [q for q in out["rows"] if not (q["config"] == name and q["rounds"] == r)] + [row]
            print(row, flush=True)
            json.dump(out, open("results/milp_rk.json", "w"), indent=1)

if __name__ == "__main__":
    rounds = [int(x) for x in sys.argv[1].split(",")]
    names = sys.argv[3].split("|") if len(sys.argv) > 3 else list(CFGS)
    run(rounds, float(sys.argv[2]), names)

def ks_only(steps, schedule, tl):
    """Min active key-schedule S-boxes over the first `steps` state updates, for a non-zero key difference."""
    g = json.load(open("results/ga_present.json")); prm = GAKSPresentParams.from_bits(g["chromosome"])
    _ZERO.clear(); m = Model(); k0 = m.var(80); act = []
    m.con({v: 1 for v in k0}, 1, np.inf)
    st = list(k0)
    for t in range(steps):
        st = gaks_step(m, st, t, prm, act) if schedule == "gaks" else present_ks_step(m, st, t + 1, act)
    t0 = time.time(); res = m.solve(tl)
    return dict(schedule=schedule, steps=steps, optimal=bool(res.status == 0), seconds=round(time.time() - t0, 1),
                best_found=float(res.fun) if res.x is not None else None,
                lower_bound=float(res.fun) if res.status == 0 else float(getattr(res, "mip_dual_bound", np.nan)))

def ks_only_des(steps, tl):
    """GAKS-DES: min active key-schedule S-boxes over the first `steps` updates, non-zero 56-bit key difference."""
    from gaks_core import final_des_params, SHIFTS
    prm = final_des_params()
    _ZERO.clear(); m = Model(); act = []
    C, D = m.var(28), m.var(28)            # C[i] = bit i of the 28-bit half
    m.con({v: 1 for v in C + D}, 1, np.inf)
    rot = lambda x, s: [x[(i - s) % 28] for i in range(28)]
    for t in range(steps):
        r = t - prm.w; s = SHIFTS[r] if r >= 0 else 1
        C, D = rot(C, s), rot(D, s)
        mk = prm.masks[t]
        for i in range(28):
            if (mk >> i) & 1: C[i], D[i] = D[i], C[i]
        D = [m.xor(D[i], C[(i - prm.alpha[t]) % 28]) for i in range(28)]
        C = [m.xor(C[i], D[(i - prm.beta[t]) % 28]) for i in range(28)]
        Sv = D + C                          # 56-bit state, bit j: j<28 -> D[j], else C[j-28]
        for j in range(14):
            if (prm.smask >> j) & 1:
                xs = Sv[4 * j:4 * j + 4]
                ys, A = m.sbox(xs); act.append(A); Sv[4 * j:4 * j + 4] = ys
        D, C = Sv[:28], Sv[28:]
    t0 = time.time(); res = m.solve(tl)
    return dict(schedule="gaks-des", steps=steps, optimal=bool(res.status == 0), seconds=round(time.time() - t0, 1),
                best_found=float(res.fun) if res.x is not None else None)
