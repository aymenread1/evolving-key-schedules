"""sensitivity.py -- fitness-weight sensitivity: evaluates the GA designs obtained with perturbed weights
(results/sensitivity/*.json) and the baseline on fresh keys; reports which fitness terms are active."""
import json, glob, os
import numpy as np
from gaks_core import *
from ga_design import bits_of, des_fitness, diffusion_metrics, CORR_KEYS, interkey_corr, NC, DES_KEYS, FLIP_DES
import ga_design as G
from linear_bypass import t_lin_des

rng = np.random.default_rng(31337)                  # fresh keys, disjoint from design/validation/evaluation seeds
N = 300
K = rng.integers(0, 2**63, N, dtype=np.uint64) * U(2)
F = np.array([U(1) << U(b) for b in EFFECTIVE_BITS], dtype=np.uint64)
ZTHR = 5.17 * 0.5 / np.sqrt(N)

def evaluate(prm):
    b = ks_gaks(K, prm); f = ks_gaks((K[:, None] ^ F[None, :]).ravel(), prm).reshape(N, 56, 16)
    x = f ^ b[:, None, :]
    a_r = (popcount(x) / 48).mean(axis=(0, 1)); d = np.abs(bits_of(x, 48).mean(0) - 0.5)
    sig = (d > ZTHR).sum(axis=(0, 2))
    return dict(cost=prm.sbox_count(), a_min=float(a_r.min()), sig_K1=int(sig[0]), sig_K2_16=int(sig[1:].sum()),
                T_lin=t_lin_des(prm), rho=interkey_corr(ks_gaks(CORR_KEYS, prm)))

def terms(prm):
    """fitness components on the design keys (which terms are active at the optimum?)"""
    base = ks_gaks(DES_KEYS, prm)
    fk = (DES_KEYS[:, None] ^ np.array(FLIP_DES, dtype=np.uint64)[None, :]).ravel()
    fl = ks_gaks(fk, prm).reshape(len(DES_KEYS), 56, 16)
    x = fl ^ base[:, None, :]
    a_r = (popcount(x) / 48).mean(axis=(0, 1)); p = bits_of(x, 48).mean(axis=0); dev = np.abs(p - 0.5)
    chi = (4 * x.shape[0] * (p - 0.5) ** 2).mean(axis=(0, 2)) - 1.0
    ic = interkey_corr(ks_gaks(CORR_KEYS, prm))
    return dict(avalanche_terms=float(np.abs(a_r - .5).mean() + np.abs(a_r - .5).max()),
                sac_terms=float(dev.mean() + dev.mean(axis=(0, 2)).max()),
                det_fraction=float((dev == 0.5).mean(axis=(0, 2)).max()),
                chi_max=float(chi.max()), chi_excess_over_0_15=float(max(0, chi.max() - 0.15)),
                rho=ic, rho_excess=float(max(0, ic - 4 / np.sqrt(NC))),
                cost_term_lambda1=prm.sbox_count() / ((GAKSParams.MAXW + 16) * 14))

out = {"baseline": {"settings": dict(lambda_=0.5, chi_tol=0.15, w_det=5, w_rho=5)}}
gd = json.load(open("results/ga_des.json")); P = GAKSParams.from_bits(gd["chromosome"])
out["baseline"]["ga_output"] = evaluate(P); out["baseline"]["active_terms"] = terms(P)
P7 = final_des_params(); out["baseline"]["final_w7"] = evaluate(P7)
for fn in sorted(glob.glob("results/sensitivity/ga_des_sens_*.json")):
    g = json.load(open(fn)); prm = GAKSParams.from_bits(g["chromosome"])
    p7 = GAKSParams(prm.masks, prm.alpha, prm.beta, prm.smask, prm.w + 1)
    name = os.path.basename(fn)[len("ga_des_sens_"):-5]
    out[name] = {"settings": dict(lambda_=g["lambda_"], chi_tol=g["chi_tol"], w_det=g["w_det"], w_rho=g["w_rho"]),
                 "w": prm.w, "sboxes_per_step": bin(prm.smask).count("1"),
                 "ga_output": evaluate(prm), "plus_one_warmup": evaluate(p7)}
    print(name, out[name]["settings"], out[name]["w"], out[name]["sboxes_per_step"], out[name]["ga_output"], out[name]["plus_one_warmup"])
print("baseline", out["baseline"])
json.dump(out, open("results/sensitivity/summary.json", "w"), indent=1)
