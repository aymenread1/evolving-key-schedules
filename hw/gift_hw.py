"""gift_hw.py -- round-based Verilog for the GIFT-64 key schedule and GAKS-GIFT (generated from the Python
parameters), bit-exact check of the GAKS-GIFT step against Python with Yosys 'eval', and synthesis with the
NanGate 45 nm library (same flow as synth_nangate.py). Writes hw/area_gift_nangate45.json."""
import json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE); os.chdir(ROOT)
import numpy as np
from gift_core import *
from gen_verilog import rotl

ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
P = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); P.w = ev["final_w"]
TOTAL = P.w + ROUNDS - 1
SBF = "  function [3:0] sbox(input [3:0] x);\n    case (x)\n" + "\n".join(f"      4'h{x:x}: sbox = 4'h{int(GS[x]):x};" for x in range(16)) + "\n    endcase\n  endfunction\n"
UPD = """  function [127:0] upd(input [127:0] k);
    upd = {k[31:16] >> 2 | k[31:16] << 14, k[15:0] >> 12 | k[15:0] << 4, k[127:32]};
  endfunction
"""

def v_gift():
    return f"""module ks_gift (input clk, input load, input [127:0] key, output [31:0] rk);
{UPD}  reg [127:0] st;
  assign rk = st[31:0];
  always @(posedge clk) st <= load ? key : upd(st);
endmodule
"""

def v_gaks():
    nib = [j for j in range(32) if (P.smask >> j) & 1]
    sb = "\n".join(f"    s[{4*j+3}:{4*j}] = sbox(s[{4*j+3}:{4*j}]);" for j in nib)
    genes = "\n".join(f"""      1'd{b}: begin
        dm = (h0 ^ l0) & 64'h{P.masks[b]:016x}; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {rotl('h1', 64, P.alpha[b])}; h2 = h1 ^ {rotl('l2', 64, P.beta[b])};
      end""" for b in range(P.nblocks))
    assert P.nblocks == 2
    return f"""module gaks_gift_step (input [127:0] si, input [5:0] t, output reg [127:0] so);
{SBF}{UPD}  reg [127:0] u, s; reg [63:0] h0, l0, dm, h1, l1, h2, l2; wire [5:0] t1 = t + 6'd1;
  always @(*) begin
    u = upd(si); h0 = u[127:64]; l0 = u[63:0]; dm = 0; h1 = 0; l1 = 0; h2 = 0; l2 = 0;
    case (t[0])
{genes}
    endcase
    s = {{h2, l2}};
{sb}
    so = s ^ {{58'd0, t1, 58'd0, t1}};
  end
endmodule

module ks_gaks_gift (input clk, input load, input [127:0] key, output [31:0] rk);
  reg [127:0] st; reg [5:0] t; wire [127:0] nx;
  gaks_gift_step u0 (.si(st), .t(t), .so(nx));
  assign rk = st[31:0];
  always @(posedge clk) begin
    if (load) begin st <= key; t <= 6'd0; end
    else if (t != 6'd{TOTAL}) begin st <= nx; t <= t + 6'd1; end
  end
endmodule
"""

def ref_step(v, t):
    H, L = as_u64(v >> 64), as_u64(v & ((1 << 64) - 1))
    H, L = gaks_gift_step(H, L, t, P)
    return (int(H[0]) << 64) | int(L[0])

if __name__ == "__main__":
    open("hw/ks_gift.v", "w").write(v_gift()); open("hw/ks_gaks_gift.v", "w").write(v_gaks())
    rng = np.random.default_rng(21)
    cases = [((int(rng.integers(0, 2**63)) << 65) ^ (int(rng.integers(0, 2**63)) << 1) ^ int(rng.integers(0, 2)), t) for t in range(TOTAL)]
    cmds = ["read_verilog hw/ks_gaks_gift.v", "synth -flatten -noabc -top gaks_gift_step"] + \
           [f"eval -set si 128'h{v:032x} -set t 6'd{t} -show so" for v, t in cases]
    out = subprocess.run([sys.executable, "-c", "import sys\nfrom pyosys import libyosys as ys\nd=ys.Design()\n[ys.run_pass(c.strip(), d) for c in sys.argv[1].split(';') if c.strip()]", "; ".join(cmds)],
                         capture_output=True, text=True).stdout
    got = [int(m, 2) for m in re.findall(r"\\so = \d+'([01]+)", out)]
    ok = len(got) == len(cases) and all(a == ref_step(v, t) for a, (v, t) in zip(got, cases))
    print("GAKS-GIFT step Verilog matches Python for all", len(cases), "steps:", ok)
    # original GIFT update: check one step via the round-key sequence of the sequential model is trivial wiring;
    # we check the update function directly
    cmds = ["read_verilog hw/_gift_upd.v", "synth -flatten -noabc -top gupd"] + [f"eval -set k 128'h{v:032x} -show o" for v, _ in cases[:8]]
    open("hw/_gift_upd.v", "w").write(f"module gupd (input [127:0] k, output [127:0] o);\n{UPD}  assign o = upd(k);\nendmodule\n")
    out = subprocess.run([sys.executable, "-c", "import sys\nfrom pyosys import libyosys as ys\nd=ys.Design()\n[ys.run_pass(c.strip(), d) for c in sys.argv[1].split(';') if c.strip()]", "; ".join(cmds)],
                         capture_output=True, text=True).stdout
    got = [int(m, 2) for m in re.findall(r"\\o = \d+'([01]+)", out)]
    def ref_upd(v):
        H, L = gift_ks_update(as_u64(v >> 64), as_u64(v & ((1 << 64) - 1))); return (int(H[0]) << 64) | int(L[0])
    ok2 = len(got) == 8 and all(a == ref_upd(v) for a, (v, _) in zip(got, cases[:8]))
    print("GIFT key update Verilog matches Python:", ok2)
    os.remove("hw/_gift_upd.v")
    # synthesis with NanGate 45 nm
    LIB = os.environ.get("NANGATE_LIB", "hw/lib/NangateOpenCellLibrary_typical.lib")
    res = {"verified_step": ok, "verified_update": ok2}
    for vf, top, lab in [("hw/ks_gift.v", "ks_gift", "GIFT-64"), ("hw/ks_gaks_gift.v", "ks_gaks_gift", "GAKS-GIFT")]:
        script = (f"read_verilog {vf}; hierarchy -top {top}; proc -norom; flatten; synth -flatten -top {top}; "
                  f"dfflibmap -liberty {LIB}; abc -liberty {LIB}; opt_clean -purge; stat -liberty {LIB}")
        out = subprocess.run([sys.executable, "-c", "import sys\nfrom pyosys import libyosys as ys\nd=ys.Design()\n[ys.run_pass(c.strip(), d) for c in sys.argv[1].split(';') if c.strip()]", script],
                             capture_output=True, text=True).stdout
        stat = out[out.rfind("Printing statistics"):]
        area = float(re.search(r"Chip area for (?:top )?module '\\?" + top + r"':\s*([\d.]+)", stat).group(1))
        cells = {m[1]: int(m[0]) for m in re.findall(r"^\s+(\d+)\s+(?:[\d.E+-]+\s+)?([A-Z][A-Z0-9_]+_X\d+)\s*$", stat, re.M)}
        ff = sum(n for c, n in cells.items() if c.startswith("DFF"))
        ffa = sum(n * {"DFF_X1": 4.522, "DFFR_X1": 5.320, "DFFS_X1": 5.320}.get(c, 0) for c, n in cells.items())
        res[lab] = {"area_um2": area, "GE": round(area / 0.798), "GE_flipflops": round(ffa / 0.798), "GE_logic": round((area - ffa) / 0.798), "flip_flops": ff}
        print(lab, res[lab])
    json.dump(res, open("hw/area_gift_nangate45.json", "w"), indent=1)
