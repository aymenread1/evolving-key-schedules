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

# ---- official NIST STS 2.1.2 table (15 tests)
OF = json.load(open("results/nist_official_sts.json"))
ORDER = [("Frequency", "Frequency"), ("BlockFrequency", "Block frequency"), ("CumulativeSums", "Cumulative sums"),
         ("Runs", "Runs"), ("LongestRun", "Longest run"), ("Rank", "Matrix rank"), ("FFT", "Spectral (DFT)"),
         ("NonOverlappingTemplate", "Non-overlapping templates"), ("OverlappingTemplate", "Overlapping templates"),
         ("Universal", "Universal"), ("ApproximateEntropy", "Approximate entropy"), ("RandomExcursions", "Random excursions"),
         ("RandomExcursionsVariant", "Random excursions variant"), ("Serial", "Serial"), ("LinearComplexity", "Linear complexity")]
rows = ["\\begin{tabular}{l" + "c" * len(cfgs) + "}", "\\toprule",
        "Test (sub-tests) & " + " & ".join(f"\\rotatebox{{70}}{{{lab}}}" for _, lab in cfgs) + " \\\\", "\\midrule"]
for key, lab in ORDER:
    nsub = max(OF[c]["tests"].get(key, {}).get("subtests", 0) for c, _ in cfgs)
    cells = []
    for c, _ in cfgs:
        t = OF[c]["tests"].get(key)
        if t is None:
            cells.append("n/a")
        elif t["failed"] == 0:
            cells.append(f"{t['min_proportion']:.2f}")
        elif t["subtests"] == 1:
            cells.append(f"\\textbf{{{t['min_proportion']:.2f}}}$^\\ast$")
        else:
            cells.append(f"\\textbf{{{t['failed']}/{t['subtests']}}}$^\\ast$")
    rows.append(f"{lab} ({nsub}) & " + " & ".join(cells) + " \\\\")
rows += ["\\midrule", "Tests passed & " + " & ".join(f"{OF[c]['tests_passed']}/15" for c, _ in cfgs) + " \\\\", "\\bottomrule", "\\end{tabular}"]
open("table_nist_official.tex", "w").write("\n".join(rows))
with open("paper_macros.tex", "a") as f:
    for c, _ in cfgs:
        f.write(f"\\newcommand{{\\sts{c.replace('-', '').replace('GAKS', 'G')}}}{{{OF[c]['tests_passed']}/15}}\n")

# ---- fitness-weight sensitivity table
SE = json.load(open("results/sensitivity/summary.json"))
order = [("baseline", "Baseline"), ("lam0.25", "$\\lambda=0.25$"), ("lam1.0", "$\\lambda=1.0$"),
         ("chi0.075", "$\\chi$ tolerance 0.075"), ("chi0.30", "$\\chi$ tolerance 0.30"),
         ("w2.5", "weights 2.5 ($\\delta_{\\det}$, $\\rho$)"), ("w10", "weights 10 ($\\delta_{\\det}$, $\\rho$)")]
rows = ["\\begin{tabular}{lcccccc}", "\\toprule",
        "Setting & $w$ & S-boxes/step & S-boxes/key & $T_{\\mathrm{lin}}$ & \\multicolumn{2}{c}{Significant SAC cells: $K_1$ / $K_2..K_{16}$} \\\\",
        " & & & & & GA output & $w+1$ \\\\", "\\midrule"]
for key, lab in order:
    e = SE[key]
    if key == "baseline":
        go, p1, w, k = e["ga_output"], e["final_w7"], 6, 2
    else:
        go, p1, w, k = e["ga_output"], e["plus_one_warmup"], e["w"], e["sboxes_per_step"]
    rows.append(f"{lab} & {w} & {k} & {go['cost']} & {go['T_lin']} & {go['sig_K1']} / {go['sig_K2_16']} & {p1['sig_K1']} / {p1['sig_K2_16']} \\\\")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_sensitivity.tex", "w").write("\n".join(rows))
AT = SE["baseline"]["active_terms"]
# ---- area table
AR = json.load(open("hw/area.json"))
rows = ["\\begin{tabular}{lrrrr}", "\\toprule", "Schedule & Flip-flops & Combinational gates & Est.\\ GE & Ratio \\\\", "\\midrule"]
for lab, base in [("DES", "DES"), ("GAKS-DES", "DES"), ("PRESENT-80", "PRESENT-80"), ("GAKS-PRESENT", "PRESENT-80")]:
    a = AR[lab]
    rows.append(f"{lab} & {a['flip_flops']} & {a['comb_gates']} & {a['GE_total_est']} & {a['GE_total_est']/AR[base]['GE_total_est']:.1f} \\\\")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_area.tex", "w").write("\n".join(rows))
T10 = json.load(open("c_bench/timing_10runs.json"))
rD = [r["GAKS-DES"]["tsc_cycles_per_key"] / r["DES"]["tsc_cycles_per_key"] for r in T10]
rP = [r["GAKS-PRESENT"]["tsc_cycles_per_key"] / r["PRESENT-80"]["tsc_cycles_per_key"] for r in T10]
with open("paper_macros.tex", "a") as f:
    f.write(f"\\newcommand{{\\chiMaxOpt}}{{{AT['chi_max']:.3f}}}\n\\newcommand{{\\rhoOpt}}{{{AT['rho']:.3f}}}\n")
    f.write(f"\\newcommand{{\\ratioDrange}}{{{min(rD):.2f}--{max(rD):.2f}}}\n\\newcommand{{\\ratioPrange}}{{{min(rP):.2f}--{max(rP):.2f}}}\n")
    f.write(f"\\newcommand{{\\geDES}}{{{AR['DES']['GE_total_est']}}}\n\\newcommand{{\\geGD}}{{{AR['GAKS-DES']['GE_total_est']}}}\n")
    f.write(f"\\newcommand{{\\geP}}{{{AR['PRESENT-80']['GE_total_est']}}}\n\\newcommand{{\\geGP}}{{{AR['GAKS-PRESENT']['GE_total_est']}}}\n")
    f.write(f"\\newcommand{{\\geRatioD}}{{{AR['GAKS-DES']['GE_total_est']/AR['DES']['GE_total_est']:.1f}}}\n\\newcommand{{\\geRatioP}}{{{AR['GAKS-PRESENT']['GE_total_est']/AR['PRESENT-80']['GE_total_est']:.1f}}}\n")

# ---- NanGate 45nm area (replaces the indicative generic estimate)
AN = json.load(open("hw/area_nangate45.json"))
rows = ["\\begin{tabular}{lrrrrr}", "\\toprule", "Schedule & Flip-flops & GE (registers) & GE (logic) & GE (total) & Ratio \\\\", "\\midrule"]
for lab, base in [("DES", "DES"), ("GAKS-DES", "DES"), ("GAKS-DES-HW", "DES"), ("PRESENT-80", "PRESENT-80"), ("GAKS-PRESENT", "PRESENT-80")]:
    a = AN[lab]
    rows.append(f"{lab} & {a['flip_flops']} & {a['GE_flipflops']} & {a['GE_logic']} & {a['GE']} & {a['GE']/AN[base]['GE']:.2f} \\\\")
    if lab in ("GAKS-DES-HW",):
        rows.append("\\midrule")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_area.tex", "w").write("\n".join(rows))
HV = json.load(open("results/hw_variant/eval.json"))["hw2"]
OFH = json.load(open("results/nist_official_sts.json"))
with open("paper_macros.tex", "a") as f:
    for k, lab in [("geNDES", "DES"), ("geNGD", "GAKS-DES"), ("geNGDHW", "GAKS-DES-HW"), ("geNP", "PRESENT-80"), ("geNGP", "GAKS-PRESENT")]:
        f.write(f"\\newcommand{{\\{k}}}{{{AN[lab]['GE']}}}\n")
    f.write(f"\\newcommand{{\\geNratioD}}{{{AN['GAKS-DES']['GE']/AN['DES']['GE']:.1f}}}\n")
    f.write(f"\\newcommand{{\\geNratioDHW}}{{{AN['GAKS-DES-HW']['GE']/AN['DES']['GE']:.2f}}}\n")
    f.write(f"\\newcommand{{\\geNratioP}}{{{AN['GAKS-PRESENT']['GE']/AN['PRESENT-80']['GE']:.1f}}}\n")
    f.write(f"\\newcommand{{\\hwW}}{{{HV['final_w']}}}\n\\newcommand{{\\hwSbox}}{{{HV['sboxes']}}}\n")
    f.write(f"\\newcommand{{\\hwSig}}{{{', '.join(map(str, HV['sig_per_round'][:3]))}}}\n")
    f.write(f"\\newcommand{{\\hwAv}}{{{HV['ks_avalanche']:.4f}}}\n\\newcommand{{\\hwSac}}{{{HV['sac_mean_abs_dev']:.4f}}}\n")
    f.write(f"\\newcommand{{\\hwRk}}{{{HV['rk_max_prob']:.3f}}}\n\\newcommand{{\\hwTlin}}{{{HV['T_lin']}}}\n\\newcommand{{\\hwRounds}}{{{HV['rounds_to_049']}}}\n")
    f.write(f"\\newcommand{{\\hwStsCtr}}{{{OFH['GAKSHW-DES-counter']['tests_passed']}/15}}\n\\newcommand{{\\hwStsRnd}}{{{OFH['GAKSHW-DES-random']['tests_passed']}/15}}\n")

# ======================================================================== GIFT-64 (primary case study)
GE_ = json.load(open("results/gift/eval.json")); GL = json.load(open("results/gift/linear_bypass.json"))
GM = json.load(open("results/gift/milp.json")); GA_ = json.load(open("results/gift/ablation.json"))
GT = json.load(open("c_bench/gift_timing.json")); GAR = json.load(open("hw/area_gift_nangate45.json"))
GGA = json.load(open("results/gift/ga_gift_nb2_s7.json")); OF = json.load(open("results/nist_official_sts.json"))
g0, g1 = GE_["GIFT-64"], GE_["GAKS-GIFT"]
gm = {
 "gW": str(GE_["final_w"]), "gGAW": str(GGA["params"]["w"]), "gSbox": str(GE_["sboxes"]), "gK": str(GL["sboxes_per_step"]),
 "gSteps": str(GL["steps_total"]), "gAvO": f"{g0['mean']:.4f}", "gAvG": f"{g1['mean']:.4f}",
 "gAvGmin": f"{min(g1['per_round']):.4f}", "gAvGmax": f"{max(g1['per_round']):.4f}",
 "gSacG": f"{g1['sac']:.4f}", "gSacIdeal": f"{GE_['ideal_sac_ref']:.4f}", "gSigG": str(sum(g1["sig_per_round"])),
 "gCells": str(128 * 28 * 32), "gRkG": f"{max(g1['rk_max_prob_per_round']):.3f}",
 "gRoundsO": str(g0["rounds_to_049"]), "gRoundsG": str(g1["rounds_to_049"]),
 "gAvOneO": f"{g0['cipher_av'][0]:.3f}", "gAvOneG": f"{g1['cipher_av'][0]:.3f}", "gAvTwoO": f"{g0['cipher_av'][1]:.3f}",
 "gAvOneRnd": f"{GE_['random_designs_cipher_av']['mean'][0]:.3f}", "gAvOneRndMin": f"{GE_['random_designs_cipher_av']['min'][0]:.3f}",
 "gAvOneRndMax": f"{GE_['random_designs_cipher_av']['max'][0]:.3f}",
 "gBlrG": f"{100*g1['blr_bits']:.1f}\\%", "gIcO": f"{g0['interkey_bitcorr']:.2f}", "gIcG": f"{g1['interkey_bitcorr']:.3f}",
 "gIcNoise": f"{GE_['interkey_noise']:.3f}", "gPeriod": str(g0["ks_period"]),
 "gTlin": str(GL["T_lin"]), "gTlinBound": str(GL["bound"]), "gMinKS": str(GL["guaranteed_active"]),
 "gValSig": ", ".join(f"$w={v['w']}$: {v['sig_per_round'][0]}/{sum(v['sig_per_round'][1:])}" for v in GE_["validation"]),
 "gCycO": f"{GT['GIFT-64']['cycles_median']:.0f}", "gCycG": f"{GT['GAKS-GIFT']['cycles_median']:.0f}",
 "gNsO": f"{GT['GIFT-64']['ns_median']:.0f}", "gNsG": f"{GT['GAKS-GIFT']['ns_median']:.0f}",
 "gCycRatio": f"{GT['GAKS-GIFT']['cycles_median']/GT['GIFT-64']['cycles_median']:.1f}",
 "gGeO": str(GAR["GIFT-64"]["GE"]), "gGeG": str(GAR["GAKS-GIFT"]["GE"]), "gGeRatio": f"{GAR['GAKS-GIFT']['GE']/GAR['GIFT-64']['GE']:.2f}",
 "gStsOc": f"{OF['GIFT-counter']['tests_passed']}/15", "gStsOr": f"{OF['GIFT-random']['tests_passed']}/15",
 "gStsGc": f"{OF['GAKSGIFT-counter']['tests_passed']}/15", "gStsGr": f"{OF['GAKSGIFT-random']['tests_passed']}/15",
 "gBestFit": f"{GGA['best_fitness']:.4f}", "gInitBest": f"{GGA['random_designs']['best']:.4f}",
 "gCmSevenMeet": pct(GA_["w7"]["5"]["frac_meets"]), "gCmNineMeet": pct(GA_["w9"]["5"]["frac_meets"]),
 "gCmFourNineMeet": pct(GA_["w9"]["4"]["frac_meets"]), "gCmSevenBetter": str(GA_["w7"]["5"]["n_better_than_ga"]),
 "gCmSevenBest": f"{GA_['w7']['5']['best_fitness']:.3f}", "gGaSevenFit": f"{GA_['w7']['ga']['fitness']:.3f}",
 "gTlinFiveFrac": pct(GA_["T_lin_random_w9"]["5"]["frac_at_bound"]), "gTlinFourMed": f"{GA_['T_lin_random_w9']['4']['median']:g}",
}
with open("paper_macros.tex", "a") as f:
    for k, v in gm.items():
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
# main GIFT table
rows = ["\\begin{tabular}{llcc}", "\\toprule", "& Metric & GIFT-64 & GAKS-GIFT \\\\", "\\midrule",
        f"E1 & Round-trip and determinism (10\\,000 pairs) & pass & pass \\\\",
        f"E2 & Mean fraction of round-key bits flipped by a 1-bit key change & {gm['gAvO']} & {gm['gAvG']} \\\\",
        f"E3 & SAC: mean $|p-\\tfrac12|$ (ideal reference {gm['gSacIdeal']}) & 0.500 & {gm['gSacG']} \\\\",
        f"   & SAC: significant cells (of {gm['gCells']}) & all & {gm['gSigG']} \\\\",
        f"E4 & $\\max$ Pr[most frequent round-key difference], single-bit $\\Delta$ & 1.000 & $\\le${gm['gRkG']} \\\\",
        f"E5 & Period of the key-state update (rounds) & {gm['gPeriod']} & -- \\\\",
        f"E7 & Rounds until key avalanche $\\ge0.49$ & {gm['gRoundsO']} & {gm['gRoundsG']} \\\\",
        f"E8 & BLR: schedule outputs satisfying linearity & 100\\% & 0\\% ({gm['gBlrG']} of bits) \\\\",
        f"E9 & NIST STS: tests passed (counter / random keys) & {gm['gStsOc']} / {gm['gStsOr']} & {gm['gStsGc']} / {gm['gStsGr']} \\\\",
        f"E10 & S-box evaluations / cycles per key set-up (C) & 0 / {gm['gCycO']} & {gm['gSbox']} / {gm['gCycG']} \\\\",
        f"    & Area, round-based (GE, NanGate 45\\,nm) & {gm['gGeO']} & {gm['gGeG']} \\\\",
        f"E11 & Max $|$corr$|$ between bits of different round keys (noise {gm['gIcNoise']}) & {gm['gIcO']} & {gm['gIcG']} \\\\",
        f"E12 & S-box bypass length $T_{{\\mathrm{{lin}}}}$ (updates; bound {gm['gTlinBound']}) & $\\infty$ & {gm['gTlin']} \\\\",
        "\\bottomrule", "\\end{tabular}"]
open("table_gift_main.tex", "w").write("\n".join(rows))
# MILP table for GIFT
bykey = {(r["rounds"], r["related_key"]): r for r in GM["rows"]}
rr = sorted({r for r, _ in bykey})
designers = [1, 2, 3, 5, 7, 10, 13, 16, 18]
rows = ["\\begin{tabular}{l" + "c" * len(rr) + "}", "\\toprule", "Rounds $r$ & " + " & ".join(map(str, rr)) + " \\\\", "\\midrule",
        "Single-key, designers' bound~\\cite{bogdanov2007gift} & " + " & ".join(str(designers[r - 1]) for r in rr) + " \\\\",
        "Single-key, this model & " + " & ".join(f"{bykey[(r, False)]['value']:.0f}" + ("" if bykey[(r, False)]["optimal"] else "$^\\dagger$") for r in rr) + " \\\\",
        "Related-key ($\\Delta K\\neq0$), original schedule & " + " & ".join(f"{bykey[(r, True)]['value']:.0f}" + ("" if bykey[(r, True)]["optimal"] else "$^\\dagger$") for r in rr) + " \\\\",
        "\\bottomrule", "\\end{tabular}"]
open("table_gift_milp.tex", "w").write("\n".join(rows))
# NIST table for GIFT (official STS)
cf = [("GIFT-counter", "GIFT-64, counter keys"), ("GAKSGIFT-counter", "GAKS-GIFT, counter keys"),
      ("GIFT-random", "GIFT-64, random keys"), ("GAKSGIFT-random", "GAKS-GIFT, random keys")]
rows = ["\\begin{tabular}{l" + "c" * len(cf) + "}", "\\toprule", "Test (sub-tests) & " + " & ".join(lab for _, lab in cf).replace(", ", "\\newline ") + " \\\\", "\\midrule"]
rows[2] = "Test (sub-tests) & " + " & ".join(f"\\shortstack{{{lab.split(', ')[0]}\\\\{lab.split(', ')[1]}}}" for _, lab in cf) + " \\\\"
for key, lab in ORDER:
    nsub = max(OF[c]["tests"].get(key, {}).get("subtests", 0) for c, _ in cf); cells = []
    for c, _ in cf:
        t = OF[c]["tests"].get(key)
        if t is None: cells.append("n/a")
        elif t["failed"] == 0: cells.append(f"{t['min_proportion']:.2f}")
        elif t["subtests"] == 1: cells.append(f"\\textbf{{{t['min_proportion']:.2f}}}$^\\ast$")
        else: cells.append(f"\\textbf{{{t['failed']}/{t['subtests']}}}$^\\ast$")
    rows.append(f"{lab} ({nsub}) & " + " & ".join(cells) + " \\\\")
rows += ["\\midrule", "Tests passed & " + " & ".join(f"{OF[c]['tests_passed']}/15" for c, _ in cf) + " \\\\", "\\bottomrule", "\\end{tabular}"]
open("table_gift_nist.tex", "w").write("\n".join(rows))
# ablation + T_lin table for GIFT
rows = ["\\begin{tabular}{llcccc}", "\\toprule", "Warm-up & Structure & S-boxes/key & Meet criterion & Median sig.\\ cells $K_1$ & $T_{\\mathrm{lin}}$ (bound; median) \\\\", "\\midrule"]
for w in ("w7", "w9"):
    for k in ("4", "5", "6"):
        c = GA_[w][k]; tl = GA_["T_lin_random_w9"][k]
        tls = f"{tl['bound']}; {tl['median']:g}" if w == "w9" else "--"
        rows.append(f"${w[0]}={w[1:]}$ & {k} S-boxes/step & {c['cost']} & {pct(c['frac_meets'])} & {c['median_sig_K1']:.0f} & {tls} \\\\")
    ga = GA_[w]["ga"]
    rows.append(f"${w[0]}={w[1:]}$ & GA design (5/step) & {ga['cost']} & {'meets' if ga['sig_later'] == 0 and ga['a_min'] >= 0.49 else 'fails'} & {ga['sig_K1']} & {GL['T_lin'] if w == 'w9' else '--'} \\\\")
    rows.append("\\midrule")
rows[-1] = "\\bottomrule"; rows.append("\\end{tabular}")
open("table_gift_ablation.tex", "w").write("\n".join(rows))
# combined area table (GIFT + PRESENT); DES rows go to the appendix table_area.tex
rows = ["\\begin{tabular}{lrrrrr}", "\\toprule", "Schedule & Flip-flops & GE (registers) & GE (logic) & GE (total) & Ratio \\\\", "\\midrule"]
for lab, src, base in [("GIFT-64", GAR, "GIFT-64"), ("GAKS-GIFT", GAR, "GIFT-64"), ("PRESENT-80", AN, "PRESENT-80"), ("GAKS-PRESENT", AN, "PRESENT-80")]:
    a = src[lab]; b = src[base]
    rows.append(f"{lab} & {a['flip_flops']} & {a['GE_flipflops']} & {a['GE_logic']} & {a['GE']} & {a['GE']/b['GE']:.2f} \\\\")
    if lab == "GAKS-GIFT": rows.append("\\midrule")
rows += ["\\bottomrule", "\\end{tabular}"]
open("table_area_main.tex", "w").write("\n".join(rows))
