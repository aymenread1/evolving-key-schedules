"""gen_rk_table.py -- LaTeX table and macros for the SAT related-key bounds (results/gift/sat_rk.json)."""
import json, math
d = json.load(open("results/gift/sat_rk.json"))
R = {(r["mode"], r["rounds"]): r for r in d["rows"]}
def cell(mode, r):
    q = R.get((mode, r))
    if q is None: return "--"
    return str(q["value"]) if q.get("optimal") else f"$\\ge${q['lower_bound']}"
designers = [1, 2, 3, 5, 7, 10, 13, 16, 18]
rows = ["\\begin{tabular}{rccccc}", "\\toprule",
        " & \\multicolumn{2}{c}{Single-key} & \\multicolumn{3}{c}{Related-key ($\\Delta K\\neq0$)} \\\\",
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
        "$r$ & designers~\\cite{bogdanov2007gift} & this model & GIFT-64 & GAKS-GIFT, from $K_1$ & GAKS-GIFT, from $K_2$ \\\\", "\\midrule"]
for r in range(1, 13):
    rows.append(f"{r} & {designers[r-1] if r <= 9 else '--'} & {cell('sk-orig', r)} & {cell('rk-orig', r)} & {cell('rk-gaks-free', r)} & {cell('rk-gaks-free@1', r)} \\\\")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_gift_rk.tex", "w").write("\n".join(rows))
lg = math.log2(16 / 6)
mac = {"rkWarm": str(R[("ks-gaks", 9)]["lower_bound"]), "rkFullBound": "43", "rkFullLog": f"{43 * lg:.1f}",
       "rkOrigFull": "26", "rkOrigFullLog": f"{26 * lg:.1f}", "rkFreeFull": "39", "rkFreeFullLog": f"{39 * lg:.1f}"}
with open("paper_macros.tex") as f: txt = f.read()
for k, v in mac.items():
    if f"\\newcommand{{\\{k}}}" not in txt: txt += f"\\newcommand{{\\{k}}}{{{v}}}\n"
open("paper_macros.tex", "w").write(txt)
print("\n".join(rows)); print(mac)
