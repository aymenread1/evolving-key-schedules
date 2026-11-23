"""synth_nangate.py -- synthesis of the five round-based key-schedule datapaths with Yosys (pyosys, with ABC)
mapped to the NanGate 45 nm Open Cell Library (typical corner). GE = area / area(NAND2_X1).
Library: NangateOpenCellLibrary_typical.lib from OpenROAD-flow-scripts (flow/platforms/nangate45/lib)."""
import json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
LIB = "lib/NangateOpenCellLibrary_typical.lib"   # not redistributed: fetch it from OpenROAD-flow-scripts,
# flow/platforms/nangate45/lib/, into hw/lib/ (SHA-256 of the copy used: see nangate_lib.sha256)
NAND2 = 0.798
TOPS = [("ks_des.v", "ks_des", "DES"), ("ks_gaks_des.v", "ks_gaks_des", "GAKS-DES"),
        ("ks_gaks_des_hw.v", "ks_gaks_des_hw", "GAKS-DES-HW"),
        ("ks_present.v", "ks_present", "PRESENT-80"), ("ks_gaks_present.v", "ks_gaks_present", "GAKS-PRESENT")]
RUN = """
import sys
from pyosys import libyosys as ys
d = ys.Design()
for c in sys.argv[1].split(";"):
    if c.strip(): ys.run_pass(c.strip(), d)
"""
def synth(vf, top):
    script = (f"read_verilog {vf}; hierarchy -top {top}; proc -norom; flatten; synth -flatten -top {top}; "
              f"dfflibmap -liberty {LIB}; abc -liberty {LIB}; opt_clean -purge; stat -liberty {LIB}")
    out = subprocess.run([sys.executable, "-c", RUN, script], capture_output=True, text=True).stdout
    stat = out[out.rfind("Printing statistics"):]
    area = float(re.search(r"Chip area for (?:top )?module '\\?" + top + r"':\s*([\d.]+)", stat).group(1))
    cells = {m[1]: int(m[0]) for m in re.findall(r"^\s+(\d+)\s+(?:[\d.E+-]+\s+)?([A-Z][A-Z0-9_]+_X\d+)\s*$", stat, re.M)}
    ff_area = sum(n * {"DFF_X1": 4.522, "DFFR_X1": 5.320, "DFFS_X1": 5.320, "DFFRS_X1": 6.384}.get(c, 0) for c, n in cells.items())
    return area, cells, ff_area, stat
res = {}
for vf, top, lab in TOPS:
    area, cells, ffa, stat = synth(vf, top)
    res[lab] = {"area_um2": area, "GE": round(area / NAND2), "GE_flipflops": round(ffa / NAND2),
                "GE_logic": round((area - ffa) / NAND2), "cells": cells,
                "flip_flops": sum(n for c, n in cells.items() if c.startswith("DFF"))}
    print(lab, res[lab]["GE"], res[lab]["GE_flipflops"], res[lab]["GE_logic"], res[lab]["flip_flops"])
res["setup"] = {"tool": "Yosys 0.71 (pyosys) with ABC", "library": "NanGate 45nm Open Cell Library, typical corner (OpenROAD-flow-scripts)",
                "GE_definition": "area / area(NAND2_X1) = area / 0.798 um^2",
                "architecture": "round-based: one state update per clock; state register, step counter, load multiplexer",
                "note": "synthesis for area only, no timing constraints, no scan; wiring-only permutations cost no area"}
json.dump(res, open("area_nangate45.json", "w"), indent=1)
