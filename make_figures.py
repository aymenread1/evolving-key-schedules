"""make_figures.py -- figures + cost/diffusion trade-off evaluation (results/summary.json)."""
import json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gaks_core import *
from ga_design import bits_of

plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 200, "savefig.bbox": "tight"})
R = json.load(open("results/eval.json"))
GD = json.load(open("results/ga_des.json")); GP = json.load(open("results/ga_present.json"))
import os; os.makedirs("figures", exist_ok=True)
C_STD, C_GA, C_RND = "#7f7f7f", "#1f4e9c", "#d08c2a"

# ---- Fig 1: GA convergence
fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
for a, g, t in [(ax[0], GD, "GAKS-DES"), (ax[1], GP, "GAKS-PRESENT")]:
    h = g["history"]
    a.plot([x["gen"] for x in h], [x["best"] for x in h], color=C_GA, label="best")
    a.plot([x["gen"] for x in h], [x["mean"] for x in h], color=C_RND, ls="--", label="population mean")
    a.axhline(g["random_designs"]["best"], color=C_STD, ls=":", label="best random initial design")
    a.set_ylim(-1.0, 0.0)
    a.set_title(t); a.set_xlabel("generation"); a.set_ylabel("fitness (clipped at -1)")
ax[0].legend(fontsize=7)
plt.savefig("figures/fig_ga_convergence.pdf"); plt.savefig("figures/fig_ga_convergence.png"); plt.close()

# ---- Fig 2: key-schedule avalanche per round
fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
r16 = np.arange(1, 17)
ax[0].plot(r16, R["E2"]["DES"]["per_round_mean"], "s-", color=C_STD, ms=3, label="DES (standard)")
ax[0].plot(r16, R["E2"]["GAKS-DES"]["per_round_mean"], "o-", color=C_GA, ms=3, label="GAKS-DES")
ax[0].set_title("DES: 1-bit key flip"); ax[0].set_xlabel("round $r$"); ax[0].set_ylabel("mean fraction of $K_r$ bits flipped")
r32 = np.arange(1, 33)
ax[1].plot(r32, R["E11"]["PRESENT-80"]["ks_per_round_mean"], "s-", color=C_STD, ms=2, label="PRESENT-80 (standard)")
ax[1].plot(r32, R["E11"]["GAKS-PRESENT"]["ks_per_round_mean"], "o-", color=C_GA, ms=2, label="GAKS-PRESENT")
ax[1].set_title("PRESENT-80: 1-bit key flip"); ax[1].set_xlabel("round $r$")
for a in ax: a.axhline(0.5, color="k", lw=0.6, ls=":"); a.set_ylim(-0.02, 0.6); a.legend(fontsize=7)
plt.savefig("figures/fig_ks_avalanche.pdf"); plt.savefig("figures/fig_ks_avalanche.png"); plt.close()

# ---- Fig 3: round-reduced cipher key avalanche
fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
ax[0].plot(r16[:8], R["E7"]["DES"]["per_round_mean"][:8], "s-", color=C_STD, ms=3, label="DES")
ax[0].plot(r16[:8], R["E7"]["GAKS-DES"]["per_round_mean"][:8], "o-", color=C_GA, ms=3, label="DES + GAKS")
rd = R["E7"]["random-designs-20"]
ax[0].fill_between(r16[:8], rd["per_round_min"][:8], rd["per_round_max"][:8], color=C_RND, alpha=0.2, lw=0)
ax[0].plot(r16[:8], rd["per_round_mean"][:8], "^--", color=C_RND, ms=3, label="DES + 20 random schedules")
ax[0].set_title("DES state after $r$ rounds"); ax[0].set_xlabel("rounds"); ax[0].set_ylabel("fraction of state bits flipped")
ax[1].plot(r32[:10], R["E11"]["PRESENT-80"]["cipher_key_avalanche_per_round"][:10], "s-", color=C_STD, ms=3, label="PRESENT-80")
ax[1].plot(r32[:10], R["E11"]["GAKS-PRESENT"]["cipher_key_avalanche_per_round"][:10], "o-", color=C_GA, ms=3, label="PRESENT-80 + GAKS")
ax[1].set_title("PRESENT state after $r$ rounds"); ax[1].set_xlabel("rounds")
for a in ax: a.axhline(0.5, color="k", lw=0.6, ls=":"); a.legend(fontsize=7)
plt.savefig("figures/fig_cipher_avalanche.pdf"); plt.savefig("figures/fig_cipher_avalanche.png"); plt.close()

# ---- Fig 4: related-key predictability
fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
ax[0].semilogy(r16, R["E4"]["DES"]["max_prob_subkey_diff_per_round"], "s-", color=C_STD, ms=3, label="DES")
ax[0].semilogy(r16, R["E4"]["GAKS-DES"]["max_prob_subkey_diff_per_round"], "o-", color=C_GA, ms=3, label="GAKS-DES")
ax[0].axhline(1 / 1000, color="k", lw=0.6, ls=":", label="resolution (1/N)")
ax[0].set_title("DES"); ax[0].set_xlabel("round $r$"); ax[0].set_ylabel(r"max$_\Delta$ Pr[$\Delta K_r$ = most frequent]")
ax[1].semilogy(r32, R["E11"]["PRESENT-80"]["rk_max_prob_per_round"], "s-", color=C_STD, ms=2, label="PRESENT-80")
ax[1].semilogy(r32, R["E11"]["GAKS-PRESENT"]["rk_max_prob_per_round"], "o-", color=C_GA, ms=2, label="GAKS-PRESENT")
ax[1].axhline(1 / 500, color="k", lw=0.6, ls=":", label="resolution (1/N)")
ax[1].set_title("PRESENT-80"); ax[1].set_xlabel("round $r$")
for a in ax: a.set_ylim(5e-4, 2); a.legend(fontsize=7)
plt.savefig("figures/fig_related_key.pdf"); plt.savefig("figures/fig_related_key.png"); plt.close()

# ---- Trade-off: cost vs diffusion for random designs and GA designs (fresh evaluation keys)
rng = np.random.default_rng(4242)
N = 300
k = rng.integers(0, 2**63, N, dtype=np.uint64) * U(2)
F = np.array([U(1) << U(b) for b in EFFECTIVE_BITS], dtype=np.uint64)
ZTHR = 5.17 * 0.5 / np.sqrt(N)       # Bonferroni (1%, 43,008 cells) per-cell threshold at N keys
def diff_of(prm, full=False):
    b = ks_gaks(k, prm); f = ks_gaks((k[:, None] ^ F[None, :]).ravel(), prm).reshape(N, 56, 16)
    x = f ^ b[:, None, :]
    a_r = (popcount(x) / 48).mean(axis=(0, 1)); d = np.abs(bits_of(x, 48).mean(0) - 0.5)
    if full:
        return float(a_r.min()), float(d.mean()), (d > ZTHR).sum(axis=(0, 2)).tolist()
    return float(a_r.min()), float(d.mean())
pts = []
for i in range(200):
    prm = GAKSParams.from_bits(np.random.default_rng(5000 + i).integers(0, 2, GAKSParams.NBITS))
    a, s = diff_of(prm); pts.append(dict(kind="random", cost=prm.sbox_count(), a_min=a, sac=s))
for fn in sorted(glob.glob("results/ga_des*.json")):
    g = json.load(open(fn)); prm = final_des_params(fn)
    a, s = diff_of(prm); pts.append(dict(kind="GA", lam=g["lambda_"], cost=prm.sbox_count(), a_min=a, sac=s))
ref = float(np.abs(rng.integers(0, 2, (N, 56, 16, 48)).mean(0) - 0.5).mean())
# cost-matched ablation: random designs sharing a warm-up count and k S-boxes per step,
# done for the raw GA output (w=6) and for the final design (w=7)
from ga_design import des_fitness
def costmatched(ref_prm, seed0):
    ga_full = diff_of(ref_prm, full=True)
    ga_fit = float(des_fitness(ref_prm)[0])
    out, kga = {}, bin(ref_prm.smask).count("1")
    for kk in (1, 2, 3):
        res = []
        for i in range(40):
            g = np.random.default_rng(seed0 + 100 * kk + i)
            bits = g.integers(0, 2, GAKSParams.NBITS)
            nib = g.choice(14, kk, replace=False); sm = sum(1 << int(j) for j in nib)
            off = (GAKSParams.MAXW + 16) * GAKSParams.STEP
            bits[off:off + 14] = [(sm >> (13 - q)) & 1 for q in range(14)]
            bits[off + 14:off + 17] = [(ref_prm.w >> (2 - q)) & 1 for q in range(3)]
            prm = GAKSParams.from_bits(bits)
            a, sc, sig = diff_of(prm, full=True)
            res.append(dict(a_min=a, sig_round1=sig[0], sig_later=int(sum(sig[1:])), fitness=float(des_fitness(prm)[0]), cost=prm.sbox_count()))
        ok = [r for r in res if r["a_min"] >= 0.49 and r["sig_later"] == 0]
        out[str(kk)] = dict(n=len(res), cost=res[0]["cost"], frac_meets_criterion=len(ok) / len(res),
                            median_sig_round1=float(np.median([r["sig_round1"] for r in res])),
                            fitness_median=float(np.median([r["fitness"] for r in res])),
                            fitness_best=float(max(r["fitness"] for r in res)),
                            n_better_than_reference=int(sum(r["fitness"] > ga_fit for r in res)) if kk == kga else None)
    ref_eval = dict(w=ref_prm.w, cost=ref_prm.sbox_count(), a_min=ga_full[0], sig_round1=ga_full[2][0],
                    sig_later=int(sum(ga_full[2][1:])), fitness=ga_fit)
    return out, ref_eval
cm6, ga6 = costmatched(GAKSParams.from_bits(GD["chromosome"]), 9000)
cm7, ga7 = costmatched(final_des_params(), 9000)
rnd = [p for p in pts if p["kind"] == "random"]
good = [p for p in rnd if p["a_min"] >= 0.49]
summary = {"cost_matched_w6": cm6, "ga_raw_w6": ga6, "cost_matched_w7": cm7, "final_w7": ga7, "tradeoff_points": pts, "ideal_sac_ref_N300": ref,
           "random_designs": {"n": len(rnd), "frac_a_min_ge_0.49": len(good) / len(rnd),
                              "cost_of_good_random_min": min(p["cost"] for p in good) if good else None,
                              "cost_of_good_random_median": float(np.median([p["cost"] for p in good])) if good else None,
                              "cost_all_median": float(np.median([p["cost"] for p in rnd]))}}
json.dump(summary, open("results/summary.json", "w"), indent=1)
fig, ax = plt.subplots(figsize=(3.5, 2.6))
ax.scatter([p["cost"] for p in rnd], [p["a_min"] for p in rnd], s=8, color=C_RND, alpha=0.6, label="random parameters (200)")
gas = {}
for p in pts:
    if p["kind"] == "GA":
        gas.setdefault(p["cost"], []).append(p)
for cost, grp in gas.items():
    ax.scatter(cost, grp[0]["a_min"], s=60, marker="*", color=C_GA, zorder=3)
    lams = ", ".join(str(g["lam"]) for g in sorted(grp, key=lambda g: g["lam"]))
    ax.annotate(f"GA design ($\\lambda$ = {lams})", (cost, grp[0]["a_min"]), textcoords="offset points", xytext=(-6, -14), fontsize=7, color=C_GA)
ax.set_xlabel("S-box evaluations per key setup"); ax.set_ylabel("worst-round key avalanche")
ax.axhline(0.5, color="k", lw=0.6, ls=":"); ax.legend(fontsize=7, loc="lower right")
plt.savefig("figures/fig_tradeoff.pdf"); plt.savefig("figures/fig_tradeoff.png"); plt.close()
print(json.dumps(summary["random_designs"], indent=1)); print(json.dumps({"w6": cm6, "ga6": ga6, "w7": cm7, "ga7": ga7}, indent=0)); print([p for p in pts if p["kind"] == "GA"], ref)
