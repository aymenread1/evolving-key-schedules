"""synth_area.py -- technology-independent synthesis of the four round-based key-schedule datapaths with
Yosys (generic gate library, no ABC: the WebAssembly build used here does not run ABC), and an indicative
gate-equivalent (GE) estimate. Writes hw/area.json.

GE factors below are ASSUMPTIONS typical of standard-cell libraries (NAND2 = 1 GE); they are not taken
from a specific foundry library and must be replaced by real synthesis numbers before publication."""
import json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
GE = {"$_NOT_": 0.67, "$_AND_": 1.33, "$_NAND_": 1.0, "$_OR_": 1.33, "$_NOR_": 1.0, "$_XOR_": 2.33,
      "$_XNOR_": 2.33, "$_ANDNOT_": 1.33, "$_ORNOT_": 1.33, "$_MUX_": 2.33}
GE_DFF = 6.0   # per flip-flop (assumption)
TOPS = [("ks_des.v", "ks_des", "DES"), ("ks_gaks_des.v", "ks_gaks_des", "GAKS-DES"),
        ("ks_present.v", "ks_present", "PRESENT-80"), ("ks_gaks_present.v", "ks_gaks_present", "GAKS-PRESENT")]

def run(vfile, top):
    script = f"read_verilog {vfile}; hierarchy -top {top}; proc -norom; flatten; synth -flatten -noabc -top {top}; opt -full; opt_clean -purge; stat"
    out = subprocess.run(["yowasp-yosys", "-p", script], capture_output=True, text=True).stdout
    stat = out[out.rfind("Printing statistics"):]
    cells = {}
    for name, n in re.findall(r"^\s+(\d+)\s+[\d.]*\s*(\$_\w+_)\s*$", stat, re.M):
        cells[n] = int(name)
    if not cells:   # older stat layout: "  $_AND_   123"
        for n, c in re.findall(r"^\s+(\$_\w+_)\s+(\d+)\s*$", stat, re.M):
            cells[n] = int(c)
    return cells, stat

res = {}
for vf, top, lab in TOPS:
    cells, stat = run(vf, top)
    ff = sum(v for k, v in cells.items() if "DFF" in k or "SDFF" in k or "DFFE" in k)
    comb = {k: v for k, v in cells.items() if k in GE}
    unknown = {k: v for k, v in cells.items() if k not in GE and not ("DFF" in k)}
    ge_comb = sum(GE[k] * v for k, v in comb.items())
    res[lab] = {"cells": cells, "flip_flops": ff, "comb_gates": sum(comb.values()),
                "GE_comb_est": round(ge_comb), "GE_ff_est": round(GE_DFF * ff), "GE_total_est": round(ge_comb + GE_DFF * ff),
                "unmapped_cells": unknown}
    print(lab, res[lab])
res["assumptions"] = {"GE_factors": GE, "GE_per_flip_flop": GE_DFF, "flow": "yowasp-yosys 0.71: proc -norom; flatten; synth -noabc; opt -full",
                      "architecture": "round-based: one state update per clock, state register + step counter + load multiplexer"}
json.dump(res, open("area.json", "w"), indent=1)
