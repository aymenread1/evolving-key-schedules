# GAKS — reproducibility package

Code, parameters and raw results behind the revised manuscript
"GAKS: A Genetic-Algorithm-Designed, Deterministic and Nonlinear Key Schedule for Block Ciphers".

## Requirements
Python 3.10+, numpy, scipy, matplotlib. (LaTeX with tikz/booktabs/natbib/cleveref to rebuild the PDF.)

## Run everything (about 25 minutes on 2 cores)
```
python3 test_vectors.py                 # DES + PRESENT-80 published test vectors
python3 nist.py                         # self-test of the SP 800-22 implementation against the spec examples
python3 ga_design.py des                # design-time GA, ~5 min (writes results/ga_des.json)
python3 ga_design.py present            # ~6 min (writes results/ga_present.json)
python3 evaluate.py                     # experiments E1-E14 -> results/eval.json
NSEQ=300 python3 run_nist.py DES-counter GAKS-DES-counter DES-random GAKS-DES-random PRESENT-counter GAKS-PRESENT-counter
python3 make_figures.py                 # figures/ + results/summary.json (cost-matched and unrestricted ablations)
python3 linear_bypass.py                # exact GF(2) S-box bypass length T_lin -> results/linear_bypass.json
python3 milp_rk.py 1,2,3,4,5,6,7 150 "PRESENT-80 single-key|PRESENT-80 related-key"   # MILP optima (needs scipy>=1.11)
(cd c_bench && python3 gen_header.py && gcc -O2 -o bench bench.c && ./bench check && ./bench time > timing.json)
python3 gen_tables.py                   # LaTeX macros/tables generated from the JSON (no hand transcription)
pdflatex paper.tex && pdflatex paper.tex   # preview build (article class)
# paper_jisa.tex = Elsevier/JISA submission version (elsarticle; compile on Overleaf). Both share
# abstract.tex, body.tex, declarations.tex, bibliography.tex, appendix.tex and the generated macros/tables.
```

## Files
| File | Content |
| --- | --- |
| `gift_core.py`, `gift_*.py` | GIFT-64-128, GAKS-GIFT, GA, evaluation, T_lin, MILP and ablation |
| `gaks_core.py` | Vectorised DES, PRESENT-80, standard key schedules, GAKS schedules |
| `ga_design.py` | Offline GA: chromosome, fitness (Eq. 7 of the paper), operators |
| `evaluate.py` | Security/diffusion experiments |
| `nist.py`, `run_nist.py` | Independent SP 800-22 tests (DFT test excluded, see paper) |
| `linear_bypass.py`, `milp_rk.py` | Exact S-box-bypass analysis and MILP models |
| `c_bench/` | C key-set-up implementations, cross-check against Python, timing |
| `results/*.json` | Raw results used in every table and figure |

Seeds: design keys 1, validation keys 777, evaluation keys 20261122, GA 2026. Final GAKS-DES = GA chromosome with w fixed to 7 (`final_des_params()`).

## Notes before submission
* Re-run the NIST battery with the official NIST `sts-2.1.2` and report all 15 tests.
* Timing figures are from vectorised Python and are indicative only; add C cycle counts.
* All random seeds are fixed in the scripts; GA results are reproducible bit-for-bit with the same numpy version.
