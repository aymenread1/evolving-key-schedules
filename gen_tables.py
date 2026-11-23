"""gen_tables.py -- writes LaTeX tables/macros straight from the JSON results (no hand transcription)."""
import json
R = json.load(open("results/eval.json")); S = json.load(open("results/summary.json"))
GD = json.load(open("results/ga_des.json")); GP = json.load(open("results/ga_present.json"))
N = lambda c: json.load(open(f"results/nist_{c}_n300.json"))
f3 = lambda x: f"{x:.3f}"; f4 = lambda x: f"{x:.4f}"; pct = lambda x: f"{100*x:.1f}\\%"

E2, E4, E5, E6, E7, E8, E10, E11 = (R[k] for k in ["E2", "E4", "E5", "E6", "E7", "E8", "E10", "E11"])
def rounds_to(v, thr=0.49):
    return next(i + 1 for i, x in enumerate(v) if x >= thr)

macros = {
 "ksAvDES": f3(E2["DES"]["overall_mean"]), "ksAvGAKS": f4(E2["GAKS-DES"]["overall_mean"]),
 "ksAvGAKSmin": f4(min(E2["GAKS-DES"]["per_round_mean"])), "ksAvGAKSmax": f4(max(E2["GAKS-DES"]["per_round_mean"])),
 "sacGAKS": f4(E2["GAKS-DES"]["sac_mean_abs_dev"]), "sacIdeal": f4(E2["ideal_random_reference"]["sac_mean_abs_dev"]),
 "sacWithinGAKS": pct(E2["GAKS-DES"]["sac_frac_within_0.05"]), "sacWithinIdeal": pct(E2["ideal_random_reference"]["sac_frac_within_0.05"]),
 "sacKoneFrac": pct(E2["GAKS-DES"]["sac_per_round_frac_gt_0.05"][0]), "sacKoneMax": f3(E2["GAKS-DES"]["sac_per_round_max"][0]),
 "sacLaterFrac": pct(max(E2["GAKS-DES"]["sac_per_round_frac_gt_0.05"][1:])),
 "rkGAKS": f3(E4["GAKS-DES"]["max_prob_over_all_rounds"]),
 "compDES": pct(E5["DES"]["holds_fraction"]), "compGAKS": pct(E5["GAKS-DES"]["holds_fraction"]),
 "avOneDES": f3(E7["DES"]["per_round_mean"][0]), "avOneGAKS": f3(E7["GAKS-DES"]["per_round_mean"][0]),
 "avOneRnd": f3(E7["random-designs-20"]["per_round_mean"][0]),
 "avOneRndMin": f3(E7["random-designs-20"]["per_round_min"][0]), "avOneRndMax": f3(E7["random-designs-20"]["per_round_max"][0]),
 "avTwoDES": f3(E7["DES"]["per_round_mean"][1]), "avTwoGAKS": f3(E7["GAKS-DES"]["per_round_mean"][1]),
 "roundsDES": str(rounds_to(E7["DES"]["per_round_mean"])), "roundsGAKS": str(rounds_to(E7["GAKS-DES"]["per_round_mean"])),
 "blrGAKSbits": pct(E8["GAKS-DES"]["mean_fraction_of_bits_satisfying"]),
 "kpsDES": f"{E10['DES_keys_per_s_vectorised']/1e6:.2f}", "kpsGAKS": f"{E10['GAKS_keys_per_s_vectorised']/1e6:.2f}",
 "kpsDrop": f"{100*(1-E10['GAKS_keys_per_s_vectorised']/E10['DES_keys_per_s_vectorised']):.0f}\\%",
 "sboxGAKS": str(E10["GAKS_sboxes"]), "stepsGAKS": str(E10["GAKS_steps"]), "wDES": str(R["params_des"]["w"]),
 "wPRESENT": str(R["params_present"]["w"]), "sboxGAKSP": str(E11["GAKS-PRESENT"]["sboxes"]),
 "pKsAv": f3(E11["PRESENT-80"]["ks_overall_mean"]), "gpKsAv": f4(E11["GAKS-PRESENT"]["ks_overall_mean"]),
 "pSac": f3(E11["PRESENT-80"]["sac_mean_abs_dev"]), "gpSac": f4(E11["GAKS-PRESENT"]["sac_mean_abs_dev"]),
 "pBlr": pct(E11["PRESENT-80"]["blr_bits_satisfying"]), "gpBlr": pct(E11["GAKS-PRESENT"]["blr_bits_satisfying"]),
 "gpRk": f3(E11["GAKS-PRESENT"]["rk_max_prob_rounds_2_32"]),
 "pAvOne": f3(E11["PRESENT-80"]["cipher_key_avalanche_per_round"][0]), "gpAvOne": f3(E11["GAKS-PRESENT"]["cipher_key_avalanche_per_round"][0]),
 "pRounds": str(rounds_to(E11["PRESENT-80"]["cipher_key_avalanche_per_round"])),
 "gpRounds": str(rounds_to(E11["GAKS-PRESENT"]["cipher_key_avalanche_per_round"])),
 "rndFracGood": pct(S["random_designs"]["frac_a_min_ge_0.49"]), "rndGoodMin": str(S["random_designs"]["cost_of_good_random_min"]),
 "rndGoodMed": f"{S['random_designs']['cost_of_good_random_median']:.0f}",
 "gaBestDES": f4(GD["best_fitness"]), "gaRndMeanDES": f3(GD["random_designs"]["mean"]),
 "gaBestP": f4(GP["best_fitness"]), "gaRndMeanP": f3(GP["random_designs"]["mean"]),
}
LB = json.load(open("results/linear_bypass.json")); CB = json.load(open("c_bench/timing.json"))
MR = json.load(open("results/milp_rk.json")); VAL = json.load(open("results/validation_w_decision.json"))
E14 = R["E14_ga_raw_w6"]
def cmrows(cm, ref):
    k = str(round(ref["cost"] / (ref["w"] + 16)))
    return k, cm[k]
k6, c6 = cmrows(S["cost_matched_w6"], S["ga_raw_w6"]); k7, c7 = cmrows(S["cost_matched_w7"], S["final_w7"])
sig = E2["GAKS-DES"]["significant_cells_per_round_bonferroni_1pct"]
macros.update({
 "kPerStep": k7, "cmN": str(c7["n"]),
 "gaRawW": str(E14["w"]), "gaRawSbox": str(E14["sboxes"]),
 "gaRawSigOne": str(E14["significant_cells_per_round_bonferroni_1pct"][0]),
 "gaRawSigTwo": str(E14["significant_cells_per_round_bonferroni_1pct"][1]),
 "valSigOneWsix": str(VAL["GA_output_w6"]["significant_cells_per_round"][0]),
 "valSigTwoWsix": str(VAL["GA_output_w6"]["significant_cells_per_round"][1]),
 "valSigOneWseven": str(VAL["w7"]["significant_cells_per_round_bonferroni_1pct"][0]),
 "cmSixMeet": pct(c6["frac_meets_criterion"]), "cmSixBetter": str(c6["n_better_than_reference"]),
 "cmSixBest": f3(c6["fitness_best"]), "gaSixFit": f3(S["ga_raw_w6"]["fitness"]),
 "cmSevenMeet": pct(c7["frac_meets_criterion"]), "cmSevenBetter": str(c7["n_better_than_reference"]),
 "sigCells": ", ".join(map(str, sig[:3])), "sigOne": str(sig[0]), "sigTwo": str(sig[1]),
 "cellsPerRound": str(E2["GAKS-DES"]["cells_per_round"]), "zBonf": f"{E2['bonferroni_z']:.2f}",
 "icDES": f"{R['E13']['DES']:.2f}", "icGAKS": f3(R["E13"]["GAKS-DES"]), "icNoise": f3(R["E13"]["noise_ref_DES_N40000"]),
 "icGP": f3(R["E13"]["GAKS-PRESENT (K1..K8)"]),
 "gaInitBestDES": f4(GD["random_designs"]["best"]), "gaInitBestP": f4(GP["random_designs"]["best"]),
 "rkGAKSroundOne": f3(E4["GAKS-DES"]["max_prob_subkey_diff_per_round"][0]),
 "tlinGD": str(LB["GAKS-DES"]["T_lin"]), "tlinGDlb": str(LB["GAKS-DES"]["lower_bound"]),
 "tlinGP": str(LB["GAKS-PRESENT"]["T_lin"]), "tlinGPlb": str(LB["GAKS-PRESENT"]["lower_bound"]),
 "tlinP": str(LB["PRESENT-80"]["T_lin"]),
 "minKSGD": str(LB["GAKS-DES"]["guaranteed_active"]), "minKSGP": str(LB["GAKS-PRESENT"]["guaranteed_active"]),
 "minKSP": str(LB["PRESENT-80"]["guaranteed_active"]),
 "tlinRangeGD": f"{min(v for v in LB['GAKS-DES']['T_lin_all_starts'] if v)}--{max(v for v in LB['GAKS-DES']['T_lin_all_starts'] if v)}",
 "cycDES": f"{CB['DES']['tsc_cycles_per_key']:.0f}", "cycGD": f"{CB['GAKS-DES']['tsc_cycles_per_key']:.0f}",
 "cycP": f"{CB['PRESENT-80']['tsc_cycles_per_key']:.0f}", "cycGP": f"{CB['GAKS-PRESENT']['tsc_cycles_per_key']:.0f}",
 "nsDES": f"{CB['DES']['ns_per_key']:.0f}", "nsGD": f"{CB['GAKS-DES']['ns_per_key']:.0f}",
 "nsP": f"{CB['PRESENT-80']['ns_per_key']:.0f}", "nsGP": f"{CB['GAKS-PRESENT']['ns_per_key']:.0f}",
 "cycRatioD": f"{CB['GAKS-DES']['tsc_cycles_per_key']/CB['DES']['tsc_cycles_per_key']:.1f}",
 "cycRatioP": f"{CB['GAKS-PRESENT']['tsc_cycles_per_key']/CB['PRESENT-80']['tsc_cycles_per_key']:.1f}",
})
# trade-off table for T_lin
tr = LB["des_tradeoff_random_designs"]
rows = ["\\begin{tabular}{rrcc}", "\\toprule", "S-boxes/step & S-boxes/key & bound $\\lceil 56/4k\\rceil$ & $T_{\\mathrm{lin}}$ (min / median of 10 random) \\\\", "\\midrule"]
for k, v in tr.items():
    rows.append(f"{k} & {v['cost']} & {v['bound']} & {v['T_lin_min']} / {v['T_lin_median']:g} \\\\")
rows += ["\\midrule", f"GA design: {LB['GAKS-DES']['sboxes_per_step']} & {R['params_des']['sboxes_per_key_setup']} & {LB['GAKS-DES']['lower_bound']} & {LB['GAKS-DES']['T_lin']} \\\\", "\\bottomrule", "\\end{tabular}"]
open("table_tlin.tex", "w").write("\n".join(rows))
# ablation table
rows = ["\\begin{tabular}{llccc}", "\\toprule", "Warm-up & Structure & S-boxes/key & Random designs meeting criterion & Median significant cells in $K_1$ \\\\", "\\midrule"]
for w, cmx, ref, lab in [(S["ga_raw_w6"]["w"], S["cost_matched_w6"], S["ga_raw_w6"], "GA output"), (S["final_w7"]["w"], S["cost_matched_w7"], S["final_w7"], "final design")]:
    for k in ("1", "2", "3"):
        c = cmx[k]
        rows.append(f"$w={w}$ & {k} S-box{'es' if k != '1' else ''}/step & {c['cost']} & {pct(c['frac_meets_criterion'])} & {c['median_sig_round1']:.0f} \\\\")
    rows.append(f"$w={w}$ & {lab} (2/step) & {ref['cost']} & {'meets' if ref['sig_later'] == 0 and ref['a_min'] >= 0.49 else 'fails'} & {ref['sig_round1']} \\\\")
    rows.append("\\midrule")
rows[-1] = "\\bottomrule"; rows.append("\\end{tabular}")
open("table_ablation.tex", "w").write("\n".join(rows))
# MILP table (PRESENT-80, proven optima)
rr = {}
for q in MR["rows"]:
    rr.setdefault(q["rounds"], {})[q["config"]] = q
rows = ["\\begin{tabular}{lccccccc}", "\\toprule", "Rounds $r$ & " + " & ".join(str(r) for r in sorted(rr)) + " \\\\", "\\midrule"]
for name, lab in [("PRESENT-80 single-key", "Single-key"), ("PRESENT-80 related-key", "Related-key ($\\Delta K\\neq0$)")]:
    rows.append(lab + " & " + " & ".join(f"{rr[r][name]['lower_bound']:.0f}" + ("" if rr[r][name]["optimal"] else "$^{\\dagger}$") for r in sorted(rr)) + " \\\\")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_milp.tex", "w").write("\n".join(rows))
with open("paper_macros.tex", "w") as f:
    for k, v in macros.items():
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

# NIST table
cfgs = [("DES-counter", "DES, counter keys"), ("GAKS-DES-counter", "GAKS-DES, counter keys"),
        ("DES-random", "DES, random keys"), ("GAKS-DES-random", "GAKS-DES, random keys"),
        ("PRESENT-counter", "PRESENT-80, counter keys"), ("GAKS-PRESENT-counter", "GAKS-PRESENT, counter keys")]
data = {c: N(c) for c, _ in cfgs}
tests = [t["test"] for t in data["DES-counter"]["tests"]]
short = {"Frequency": "Freq.", "BlockFrequency": "Block freq.", "CumulativeSums(fwd)": "CuSum (f)",
         "CumulativeSums(rev)": "CuSum (r)", "Runs": "Runs", "LongestRun": "Longest run",
         "ApproximateEntropy": "ApEn", "Serial-1": "Serial 1", "Serial-2": "Serial 2"}
lines = ["\\begin{tabular}{l" + "c" * len(cfgs) + "}", "\\toprule",
         "Test & " + " & ".join(f"\\rotatebox{{70}}{{{lab}}}" for _, lab in cfgs) + " \\\\", "\\midrule"]
for i, t in enumerate(tests):
    cells = []
    for c, _ in cfgs:
        r = data[c]["tests"][i]
        s = f"{r['proportion']:.2f}"
        cells.append(s if r["passed"] else f"\\textbf{{{s}}}$^\\ast$")
    lines.append(f"{short[t]} & " + " & ".join(cells) + " \\\\")
lines += ["\\midrule", "Passed & " + " & ".join(f"{data[c]['passed_count']}/{len(tests)}" for c, _ in cfgs) + " \\\\",
          "\\bottomrule", "\\end{tabular}"]
open("table_nist.tex", "w").write("\n".join(lines))
with open("paper_macros.tex", "a") as f:
    f.write(f"\\newcommand{{\\nistPmin}}{{{data['DES-counter']['min_proportion']:.4f}}}\n")
    f.write(f"\\newcommand{{\\nistN}}{{{data['DES-counter']['nseq']}}}\n")

# parameter appendix
pd, pp = R["params_des"], R["params_present"]
rows = ["\\begin{tabular}{rccc}", "\\toprule", "Step $t$ & $M_t$ (hex) & $\\alpha_t$ & $\\beta_t$ \\\\", "\\midrule"]
for t, (m, a, b) in enumerate(zip(pd["masks"], pd["alpha"], pd["beta"])):
    tag = f"{t+1}" + (" (warm-up)" if t < pd["w"] else f" ($K_{{{t - pd['w'] + 1}}}$)")
    rows.append(f"{tag} & \\texttt{{{m}}} & {a} & {b} \\\\")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_params_des.tex", "w").write("\n".join(rows))
with open("paper_macros.tex", "a") as f:
    f.write(f"\\newcommand{{\\smaskDES}}{{\\texttt{{0x{pd['smask']}}}}}\n")
    f.write(f"\\newcommand{{\\smaskP}}{{\\texttt{{0x{pp['smask']}}}}}\n")
    f.write(f"\\newcommand{{\\masksP}}{{\\texttt{{{', '.join(pp['masks'])}}}}}\n")
    f.write(f"\\newcommand{{\\alphaP}}{{{', '.join(map(str, pp['alpha']))}}}\n")
    f.write(f"\\newcommand{{\\betaP}}{{{', '.join(map(str, pp['beta']))}}}\n")
print(open("paper_macros.tex").read()); print(open("table_nist.tex").read())
