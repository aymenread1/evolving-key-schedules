"""
gen_verilog.py -- generates round-based (one state update per clock) Verilog for the four key schedules,
from the same parameters as the Python reference (no hand transcription), plus combinational step modules
used to cross-check the Verilog against Python with yosys 'eval'.
Bit numbering follows gaks_core.py: bit 0 = least significant bit.
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); os.chdir(ROOT)
from gaks_core import *

P = final_des_params()
GP = GAKSPresentParams.from_bits(json.load(open("results/ga_present.json"))["chromosome"])
SB = [int(x) for x in PRESENT_S]

def sbox_func():
    cases = "\n".join(f"      4'h{x:x}: sbox = 4'h{SB[x]:x};" for x in range(16))
    return f"""  function [3:0] sbox(input [3:0] x);
    case (x)
{cases}
    endcase
  endfunction
"""

def perm_assign(dst, src, table, nin, nout):
    # FIPS numbering: output bit (nout-1-k) = input bit (nin - table[k])
    return "\n".join(f"  assign {dst}[{nout-1-k}] = {src}[{nin - p}];" for k, p in enumerate(table))

def rotl(expr, w, s):
    s %= w
    return expr if s == 0 else f"{{{expr}[{w-1-s}:0], {expr}[{w-1}:{w-s}]}}"

# ------------------------------------------------------------------ DES original
def des_orig():
    sh = "\n".join(f"      4'd{r}: two = 1'b{1 if SHIFTS[r] == 2 else 0};" for r in range(16))
    return f"""module ks_des (input clk, input load, input [63:0] key, output [47:0] rk);
  reg [27:0] C, D; reg [3:0] r; reg two;
  wire [55:0] k56; wire [55:0] st = {{C, D}};
{perm_assign('k56', 'key', PC1, 64, 56)}
{perm_assign('rk', 'st', PC2, 56, 48)}
  always @(*) begin
    case (r)
{sh}
      default: two = 1'b0;
    endcase
  end
  wire [27:0] C1 = two ? {rotl('C', 28, 2)} : {rotl('C', 28, 1)};
  wire [27:0] D1 = two ? {rotl('D', 28, 2)} : {rotl('D', 28, 1)};
  always @(posedge clk) begin
    if (load) begin C <= k56[55:28]; D <= k56[27:0]; r <= 4'd0; end
    else begin C <= C1; D <= D1; r <= r + 4'd1; end
  end
endmodule
"""

# ------------------------------------------------------------------ GAKS-DES
def gaks_des_step_body(nsteps):
    lines = []
    for t in range(nsteps):
        r = t - P.w; s = SHIFTS[r] if r >= 0 else 1
        lines.append(f"""      5'd{t}: begin
        c1 = {rotl('ci', 28, s)}; d1 = {rotl('di', 28, s)};
        dm = (c1 ^ d1) & 28'h{P.masks[t]:07x}; c2 = c1 ^ dm; d2 = d1 ^ dm;
        d3 = d2 ^ {rotl('c2', 28, P.alpha[t])}; c3 = c2 ^ {rotl('d3', 28, P.beta[t])};
        rc = 56'h{RC_DES[t]:014x};
      end""")
    return "\n".join(lines)

def gaks_des_step_module():
    n = P.w + 16
    nib = [j for j in range(14) if (P.smask >> j) & 1]
    sbl = "\n".join(f"    s[{4*j+3}:{4*j}] = sbox(s[{4*j+3}:{4*j}]);" for j in nib)
    return f"""module gaks_des_step (input [55:0] si, input [4:0] t, output reg [55:0] so);
{sbox_func()}
  reg [27:0] ci, di, c1, d1, dm, c2, d2, c3, d3; reg [55:0] rc, s;
  always @(*) begin
    ci = si[55:28]; di = si[27:0];
    c1 = 0; d1 = 0; dm = 0; c2 = 0; d2 = 0; c3 = 0; d3 = 0; rc = 0;
    case (t)
{gaks_des_step_body(n)}
      default: begin c3 = ci; d3 = di; rc = 0; end
    endcase
    s = {{c3, d3}};
{sbl}
    so = s ^ rc;
  end
endmodule
"""

def gaks_des_top():
    n = P.w + 16
    return f"""module ks_gaks_des (input clk, input load, input [63:0] key, output [47:0] rk);
  reg [55:0] st; reg [4:0] t; wire [55:0] k56, nx;
{perm_assign('k56', 'key', PC1, 64, 56)}
{perm_assign('rk', 'st', PC2, 56, 48)}
  gaks_des_step u (.si(st), .t(t), .so(nx));
  always @(posedge clk) begin
    if (load) begin st <= k56; t <= 5'd0; end
    else if (t != 5'd{n}) begin st <= nx; t <= t + 5'd1; end
  end
endmodule
"""

# ------------------------------------------------------------------ PRESENT-80 original
def present_orig():
    return f"""module ks_present (input clk, input load, input [79:0] key, output [63:0] rk);
{sbox_func()}
  reg [79:0] st; reg [4:0] i;
  wire [79:0] r61 = {rotl('st', 80, 61)};
  wire [79:0] s1 = {{sbox(r61[79:76]), r61[75:0]}};
  wire [79:0] nx = {{s1[79:20], s1[19:15] ^ (i + 5'd1), s1[14:0]}};
  assign rk = st[79:16];
  always @(posedge clk) begin
    if (load) begin st <= key; i <= 5'd0; end
    else begin st <= nx; i <= i + 5'd1; end
  end
endmodule
"""

# ------------------------------------------------------------------ GAKS-PRESENT
def gaks_present_step_module():
    total = GP.w + 31
    genes = "\n".join(f"""      2'd{g}: begin
        dm = (h0 ^ l0) & 40'h{GP.masks[g]:010x}; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {rotl('h1', 40, GP.alpha[g])}; h2 = h1 ^ {rotl('l2', 40, GP.beta[g])};
      end""" for g in range(4))
    rcs = "\n".join(f"      6'd{t}: rc = 80'h{RC_PRESENT[t]:020x};" for t in range(total))
    nibH = [j for j in range(10) if (GP.smask >> (10 + j)) & 1]; nibL = [j for j in range(10) if (GP.smask >> j) & 1]
    sbl = "\n".join([f"    h2[{4*j+3}:{4*j}] = sbox(h2[{4*j+3}:{4*j}]);" for j in nibH] +
                    [f"    l2[{4*j+3}:{4*j}] = sbox(l2[{4*j+3}:{4*j}]);" for j in nibL])
    return f"""module gaks_present_step (input [79:0] si, input [5:0] t, output reg [79:0] so);
{sbox_func()}
  reg [79:0] r61, rc; reg [39:0] h0, l0, dm, h1, l1, h2, l2;
  always @(*) begin
    r61 = {rotl('si', 80, 61)};
    h0 = r61[79:40]; l0 = r61[39:0];
    dm = 0; h1 = 0; l1 = 0; h2 = 0; l2 = 0;
    case (t[1:0])
{genes}
    endcase
{sbl}
    case (t)
{rcs}
      default: rc = 80'h0;
    endcase
    so = {{h2, l2}} ^ rc;
  end
endmodule
"""

def gaks_present_top():
    total = GP.w + 31
    return f"""module ks_gaks_present (input clk, input load, input [79:0] key, output [63:0] rk);
  reg [79:0] st; reg [5:0] t; wire [79:0] nx;
  gaks_present_step u (.si(st), .t(t), .so(nx));
  assign rk = st[79:16];
  always @(posedge clk) begin
    if (load) begin st <= key; t <= 6'd0; end
    else if (t != 6'd{total}) begin st <= nx; t <= t + 6'd1; end
  end
endmodule
"""

if __name__ == "__main__":
    os.makedirs("hw", exist_ok=True)
    open("hw/ks_des.v", "w").write(des_orig())
    open("hw/ks_gaks_des.v", "w").write(gaks_des_step_module() + "\n" + gaks_des_top())
    open("hw/ks_present.v", "w").write(present_orig())
    open("hw/ks_gaks_present.v", "w").write(gaks_present_step_module() + "\n" + gaks_present_top())
    print("written", P.w, GP.w)

# ------------------------------------------------------------------ hardware-friendly GAKS-DES (cyclic gene blocks, counter constants)
def load_hw(name="hw2"):
    g = json.load(open(f"results/hw_variant/ga_des_{name}.json"))
    prm = GAKSParams.from_bits(g["chromosome"], nblocks=g["nblocks"], rc=RC_CTR_DES)
    ev = json.load(open("results/hw_variant/eval.json"))[name]
    prm.w = ev["final_w"]
    return prm, g["nblocks"]

def gaks_des_hw(name="hw2"):
    prm, nb = load_hw(name)
    n = prm.w + 16
    bits = max(1, (nb - 1).bit_length())
    shtab = "\n".join(f"      5'd{t}: two = 1'b{1 if (t - prm.w >= 0 and SHIFTS[t - prm.w] == 2) else 0};" for t in range(n))
    genes = "\n".join(f"""      {bits}'d{b}: begin
        dm = (c1 ^ d1) & 28'h{prm.masks[b]:07x}; c2 = c1 ^ dm; d2 = d1 ^ dm;
        d3 = d2 ^ {rotl('c2', 28, prm.alpha[b])}; c3 = c2 ^ {rotl('d3', 28, prm.beta[b])};
      end""" for b in range(nb))
    nib = [j for j in range(14) if (prm.smask >> j) & 1]
    sbl = "\n".join(f"    s[{4*j+3}:{4*j}] = sbox(s[{4*j+3}:{4*j}]);" for j in nib)
    # block index = t mod nb; nb is a power of two here
    assert nb & (nb - 1) == 0
    return f"""module gaks_des_hw_step (input [55:0] si, input [4:0] t, output reg [55:0] so);
{sbox_func()}
  reg [27:0] ci, di, c1, d1, dm, c2, d2, c3, d3; reg [55:0] s; reg two; wire [4:0] t1 = t + 5'd1;
  always @(*) begin
    case (t)
{shtab}
      default: two = 1'b0;
    endcase
    ci = si[55:28]; di = si[27:0];
    c1 = two ? {rotl('ci', 28, 2)} : {rotl('ci', 28, 1)};
    d1 = two ? {rotl('di', 28, 2)} : {rotl('di', 28, 1)};
    dm = 0; c2 = 0; d2 = 0; c3 = 0; d3 = 0;
    case (t[{bits - 1}:0])
{genes}
    endcase
    s = {{c3, d3}};
{sbl}
    so = s ^ {{23'd0, t1, 23'd0, t1}};
  end
endmodule

module ks_gaks_des_hw (input clk, input load, input [63:0] key, output [47:0] rk);
  reg [55:0] st; reg [4:0] t; wire [55:0] k56, nx;
{perm_assign('k56', 'key', PC1, 64, 56)}
{perm_assign('rk', 'st', PC2, 56, 48)}
  gaks_des_hw_step u (.si(st), .t(t), .so(nx));
  always @(posedge clk) begin
    if (load) begin st <= k56; t <= 5'd0; end
    else if (t != 5'd{n}) begin st <= nx; t <= t + 5'd1; end
  end
endmodule
"""

if __name__ == "__main__" and os.path.exists("results/hw_variant/eval.json"):
    open("hw/ks_gaks_des_hw.v", "w").write(gaks_des_hw("hw2"))
    print("written hw variant")
