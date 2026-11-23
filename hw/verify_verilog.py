"""verify_verilog.py -- checks the generated Verilog against the Python reference with yosys 'eval'.
Checks: (1) GAKS-DES step for every t, (2) GAKS-PRESENT step for every t, (3) DES K1 via PC-1/rotate/PC-2 wiring,
(4) PRESENT-80 original state update. Random 56/80-bit inputs, fixed seed."""
import os, re, subprocess, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE); os.chdir(ROOT)
import numpy as np
from gaks_core import *
import gen_verilog as GV

P, GP = GV.P, GV.GP
rng = np.random.default_rng(5)
M28i = (1 << 28) - 1; M40i = (1 << 40) - 1
rot = lambda x, s, w: ((x << s) | (x >> (w - s))) & ((1 << w) - 1)
SBi = [int(v) for v in PRESENT_S]
def sbl(x, nibs):
    for j in nibs:
        x = (x & ~(0xF << 4 * j)) | (SBi[(x >> 4 * j) & 0xF] << 4 * j)
    return x

def ref_gaks_des(S, t):
    C, D = S >> 28, S & M28i
    r = t - P.w; s = SHIFTS[r] if r >= 0 else 1
    C, D = rot(C, s, 28), rot(D, s, 28)
    d = (C ^ D) & P.masks[t]; C ^= d; D ^= d
    D ^= rot(C, P.alpha[t], 28); C ^= rot(D, P.beta[t], 28)
    S = (C << 28) | D
    S = sbl(S, [j for j in range(14) if (P.smask >> j) & 1])
    return S ^ RC_DES[t]

def ref_gaks_present(S, t):
    V = rot(S, 61, 80); H, L = V >> 40, V & M40i
    g = t % 4; d = (H ^ L) & GP.masks[g]; H ^= d; L ^= d
    L ^= rot(H, GP.alpha[g], 40); H ^= rot(L, GP.beta[g], 40)
    H = sbl(H, [j for j in range(10) if (GP.smask >> (10 + j)) & 1]); L = sbl(L, [j for j in range(10) if (GP.smask >> j) & 1])
    return ((H << 40) | L) ^ RC_PRESENT[t]

def yosys_eval(vfile, top, sets, show):
    cmds = [f"read_verilog {vfile}", f"synth -flatten -noabc -top {top}"]
    for st in sets:
        cmds.append("eval " + " ".join(f"-set {k} {v}" for k, v in st) + f" -show {show}")
    out = subprocess.run(["yowasp-yosys", "-p", "; ".join(cmds)], capture_output=True, text=True).stdout
    return [int(m.replace("_", ""), 2) for m in re.findall(r"\\" + show + r" = \d+'([01_]+)", out)]

def check(name, vfile, top, cases, show, ref):
    sets = [[(k, v) for k, v in c[0]] for c in cases]
    got = yosys_eval(vfile, top, sets, show)
    exp = [ref(*c[1]) for c in cases]
    ok = len(got) == len(exp) and all(a == b for a, b in zip(got, exp))
    print(f"{name}: {len(got)}/{len(exp)} evaluated, {'MATCH' if ok else 'MISMATCH'}")
    return ok

if __name__ == "__main__":
    res = {}
    cases = []
    for t in range(P.w + 16):
        S = int(rng.integers(0, 2**56, dtype=np.uint64))
        cases.append(([("si", f"56'h{S:014x}"), ("t", f"5'd{t}")], (S, t)))
    res["gaks_des_step"] = check("GAKS-DES step (all t)", "hw/ks_gaks_des.v", "gaks_des_step", cases, "so", ref_gaks_des)
    cases = []
    for t in range(GP.w + 31):
        S = int(rng.integers(0, 2**63, dtype=np.uint64)) << 17 | int(rng.integers(0, 2**17))
        cases.append(([("si", f"80'h{S:020x}"), ("t", f"6'd{t}")], (S, t)))
    res["gaks_present_step"] = check("GAKS-PRESENT step (all t)", "hw/ks_gaks_present.v", "gaks_present_step", cases, "so", ref_gaks_present)
    # DES K1 wiring check and PRESENT-80 original update, via small wrapper modules
    open("hw/_verify_wrap.v", "w").write(f"""module des_k1 (input [63:0] key, output [47:0] rk);
  wire [55:0] k56; wire [55:0] st;
{GV.perm_assign('k56', 'key', PC1, 64, 56)}
  assign st = {{{GV.rotl('k56[55:28]', 28, 1).replace('k56[55:28][', 'k56[').replace(':0]', ':28]') if False else ''}}};
endmodule
""")
    # simpler: build the wrapper explicitly
    open("hw/_verify_wrap.v", "w").write(f"""module des_k1 (input [63:0] key, output [47:0] rk);
  wire [55:0] k56; wire [27:0] C = k56[55:28], D = k56[27:0];
  wire [55:0] st = {{{GV.rotl('C', 28, 1)}, {GV.rotl('D', 28, 1)}}};
{GV.perm_assign('k56', 'key', PC1, 64, 56)}
{GV.perm_assign('rk', 'st', PC2, 56, 48)}
endmodule
module present_step (input [79:0] st, input [4:0] i, output [79:0] nx);
{GV.sbox_func()}
  wire [79:0] r61 = {GV.rotl('st', 80, 61)};
  wire [79:0] s1 = {{sbox(r61[79:76]), r61[75:0]}};
  assign nx = {{s1[79:20], s1[19:15] ^ (i + 5'd1), s1[14:0]}};
endmodule
""")
    keys = [int(rng.integers(0, 2**63, dtype=np.uint64)) * 2 + 1 for _ in range(8)]
    res["des_k1"] = check("DES K1 (PC-1, rotation, PC-2 wiring)", "hw/_verify_wrap.v", "des_k1",
                          [([("key", f"64'h{k:016x}")], (k,)) for k in keys], "rk", lambda k: int(ks_des(as_u64(k))[0, 0]))
    def ref_present(S, i):
        V = rot(S, 61, 80); V = (V & ~(0xF << 76)) | (SBi[(V >> 76) & 0xF] << 76)
        return V ^ ((i + 1) << 15)
    cases = []
    for i in range(8):
        S = int(rng.integers(0, 2**63, dtype=np.uint64)) << 17 | int(rng.integers(0, 2**17))
        cases.append(([("st", f"80'h{S:020x}"), ("i", f"5'd{i}")], (S, i)))
    res["present_step"] = check("PRESENT-80 original update", "hw/_verify_wrap.v", "present_step", cases, "nx", ref_present)
    json.dump(res, open("hw/verification.json", "w"), indent=1)
