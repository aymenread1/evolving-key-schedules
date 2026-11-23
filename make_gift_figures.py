"""make_gift_figures.py -- figures for the GIFT-64 / PRESENT-80 version of the article."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "savefig.bbox": "tight"})
C_STD, C_GA, C_RND = "#7f7f7f", "#1f4e9c", "#d08c2a"
G = json.load(open("results/gift/eval.json")); R = json.load(open("results/eval.json"))
GG = json.load(open("results/gift/ga_gift_nb2_s7.json")); GP = json.load(open("results/ga_present.json"))

fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
for a, g, t in [(ax[0], GG, "GAKS-GIFT"), (ax[1], GP, "GAKS-PRESENT")]:
    h = g["history"]
    a.plot([x["gen"] for x in h], [x["best"] for x in h], color=C_GA, label="best")
    a.plot([x["gen"] for x in h], [x["mean"] for x in h], color=C_RND, ls="--", label="population mean")
    a.axhline(g["random_designs"]["best"], color=C_STD, ls=":", label="best random initial design")
    a.set_ylim(-1.0, 0.0); a.set_title(t); a.set_xlabel("generation"); a.set_ylabel("fitness (clipped at -1)")
ax[0].legend(fontsize=7)
plt.savefig("figures/fig_ga_convergence_gift.pdf"); plt.savefig("figures/fig_ga_convergence_gift.png"); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
r28 = np.arange(1, 29)
ax[0].plot(r28, G["GIFT-64"]["per_round"], "s-", color=C_STD, ms=2, label="GIFT-64 (original)")
ax[0].plot(r28, G["GAKS-GIFT"]["per_round"], "o-", color=C_GA, ms=2, label="GAKS-GIFT")
ax[0].set_title("GIFT-64: 1-bit key flip"); ax[0].set_xlabel("round $r$"); ax[0].set_ylabel("mean fraction of $K_r$ bits flipped")
r32 = np.arange(1, 33)
ax[1].plot(r32, R["E11"]["PRESENT-80"]["ks_per_round_mean"], "s-", color=C_STD, ms=2, label="PRESENT-80 (original)")
ax[1].plot(r32, R["E11"]["GAKS-PRESENT"]["ks_per_round_mean"], "o-", color=C_GA, ms=2, label="GAKS-PRESENT")
ax[1].set_title("PRESENT-80: 1-bit key flip"); ax[1].set_xlabel("round $r$")
for a in ax: a.axhline(0.5, color="k", lw=0.6, ls=":"); a.set_ylim(-0.02, 0.6); a.legend(fontsize=7)
plt.savefig("figures/fig_ks_avalanche_gift.pdf"); plt.savefig("figures/fig_ks_avalanche_gift.png"); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
n = 12; rr = np.arange(1, n + 1); rd = G["random_designs_cipher_av"]
ax[0].fill_between(rr, rd["min"][:n], rd["max"][:n], color=C_RND, alpha=0.2, lw=0)
ax[0].plot(rr, rd["mean"][:n], "^--", color=C_RND, ms=3, label="GIFT-64 + 10 random schedules")
ax[0].plot(rr, G["GIFT-64"]["cipher_av"][:n], "s-", color=C_STD, ms=3, label="GIFT-64")
ax[0].plot(rr, G["GAKS-GIFT"]["cipher_av"][:n], "o-", color=C_GA, ms=3, label="GIFT-64 + GAKS")
ax[0].set_title("GIFT-64 state after $r$ rounds"); ax[0].set_xlabel("rounds"); ax[0].set_ylabel("fraction of state bits flipped")
ax[1].plot(r32[:10], R["E11"]["PRESENT-80"]["cipher_key_avalanche_per_round"][:10], "s-", color=C_STD, ms=3, label="PRESENT-80")
ax[1].plot(r32[:10], R["E11"]["GAKS-PRESENT"]["cipher_key_avalanche_per_round"][:10], "o-", color=C_GA, ms=3, label="PRESENT-80 + GAKS")
ax[1].set_title("PRESENT state after $r$ rounds"); ax[1].set_xlabel("rounds")
for a in ax: a.axhline(0.5, color="k", lw=0.6, ls=":"); a.legend(fontsize=7)
plt.savefig("figures/fig_cipher_avalanche_gift.pdf"); plt.savefig("figures/fig_cipher_avalanche_gift.png"); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(7, 2.4))
ax[0].semilogy(r28, G["GIFT-64"]["rk_max_prob_per_round"], "s-", color=C_STD, ms=2, label="GIFT-64")
ax[0].semilogy(r28, G["GAKS-GIFT"]["rk_max_prob_per_round"], "o-", color=C_GA, ms=2, label="GAKS-GIFT")
ax[0].axhline(1 / 500, color="k", lw=0.6, ls=":", label="resolution (1/N)"); ax[0].set_title("GIFT-64")
ax[0].set_xlabel("round $r$"); ax[0].set_ylabel(r"max$_\Delta$ Pr[$\Delta K_r$ = most frequent]")
ax[1].semilogy(r32, R["E11"]["PRESENT-80"]["rk_max_prob_per_round"], "s-", color=C_STD, ms=2, label="PRESENT-80")
ax[1].semilogy(r32, R["E11"]["GAKS-PRESENT"]["rk_max_prob_per_round"], "o-", color=C_GA, ms=2, label="GAKS-PRESENT")
ax[1].axhline(1 / 500, color="k", lw=0.6, ls=":", label="resolution (1/N)"); ax[1].set_title("PRESENT-80"); ax[1].set_xlabel("round $r$")
for a in ax: a.set_ylim(5e-4, 2); a.legend(fontsize=7)
plt.savefig("figures/fig_related_key_gift.pdf"); plt.savefig("figures/fig_related_key_gift.png"); plt.close()
print("ok")
