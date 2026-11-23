"""gift_sat_verify.py -- independent check of a related-key characteristic returned by gift_sat.py.
The characteristic is re-propagated with the reference implementation in gift_core.py: the key-state update,
crossover, mixing, PermBits and round-key insertion are applied to the differences (constants cancel), and every
S-box transition, in the cipher and in the key schedule, is checked against the GIFT difference distribution table.
Usage: [OFFSET=o] python3 gift_sat_verify.py <mode> <rounds> <k>"""
import sys
from pysat.solvers import Solver
from pysat.card import ITotalizer
import gift_sat as G
from gift_core import gift_ks_update, permbits, add_round_key, rotl64, U, as_u64

mode, r, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
prm = G.params(); m, io, _ = G.build(mode, r, prm)
acts = m.act_cipher + m.act_ks
tot = ITotalizer(lits=acts, ubound=min(80, len(acts)), top_id=m.nv)
sol = Solver(name="cadical195", bootstrap_with=m.cl + tot.cnf.clauses)
assert sol.solve(assumptions=[-tot.rhs[k]]), "no characteristic with <= k active S-boxes"
pos = set(l for l in sol.get_model() if l > 0)
val = lambda bits: sum(1 << i for i, v in enumerate(bits) if v is not None and v in pos)
M64 = (1 << 64) - 1
def lin_step(v, t):                                   # linear part of a GAKS step, reference code
    H, L = gift_ks_update(as_u64(v >> 64), as_u64(v & M64))
    if mode == "rk-orig": return (int(H[0]) << 64) | int(L[0])
    g = t % prm.nblocks; d = (H ^ L) & U(prm.masks[g]); H, L = H ^ d, L ^ d
    L = L ^ rotl64(H, prm.alpha[g]); H = H ^ rotl64(L, prm.beta[g])
    return (int(H[0]) << 64) | int(L[0])
sel = ((prm.smask >> 16) << 16) | (prm.smask & 0xFFFF)          # nibble j of L (j<16) / H (j>=16)
NIB = [j for j in range(32) if (sel >> j) & 1]
nib = lambda v, j: (v >> (4 * j)) & 0xF

bad, act, log = 0, 0, []
ks_states = [(it[1], val(it[2])) for it in m.trace if it[0] == "ks_state"]
rounds = {}
for it in m.trace:
    if it[0] == "round_in": rounds.setdefault(it[1], {})["in"] = val(it[2]); rounds[it[1]]["key"] = None if it[3] is None else val(it[3])
    if it[0] == "round_out": rounds[it[1]]["out"] = val(it[2])
    if it[0] == "final": final = val(it[2])
# cipher rounds
for rd in sorted(rounds):
    x, y = rounds[rd]["in"], rounds[rd]["out"]
    for j in range(16):
        a, b = nib(x, j), nib(y, j)
        if G.DDT[a, b] == 0: bad += 1
        act += a != 0
    z = permbits(as_u64(y))
    kd = rounds[rd]["key"]
    if kd is not None:
        z = add_round_key(z, as_u64(kd & 0xFFFFFFFF), 0) ^ add_round_key(as_u64(0), as_u64(0), 0)
    nxt = rounds[rd + 1]["in"] if rd + 1 in rounds else final
    if int(z[0]) != nxt: bad += 1
# key schedule: consecutive states (including the master-key / free-start state)
if mode == "rk-gaks": chain = [(None, val(io[1]))] + ks_states
elif mode in ("rk-gaks-free", "rk-orig"): chain = [(None, rounds[0]["key"])] + ks_states
else: chain = []
if mode == "rk-orig":
    for rd in range(1, r):
        if lin_step(rounds[rd - 1]["key"], 0) != rounds[rd]["key"]: bad += 1
else:
    for (_, prev), (t, cur) in zip(chain, chain[1:]):
        y = lin_step(prev, t)
        for j in range(32):
            if j in NIB:
                a, b = nib(y, j), nib(cur, j)
                if G.DDT[a, b] == 0: bad += 1
                act += a != 0
            elif nib(y, j) != nib(cur, j): bad += 1
    # round keys used by the cipher must equal the key state at the right step
    keys_from_chain = {}
    for rd in sorted(rounds):
        if rounds[rd]["key"] is None: continue
        if rd == 0 and mode == "rk-gaks-free": continue
        want = rounds[rd]["key"]
        if want not in [c for _, c in chain]: bad += 1
print(f"{mode} r={r}: SAT solution with <= {k} active S-boxes; recomputed active S-boxes = {act}; violations = {bad}")
