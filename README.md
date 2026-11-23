# Evolving Key Schedules Offline

Code, parameters, seeds and raw results for the article
**"Evolving Key Schedules Offline: A Cost-Aware Genetic Design of Deterministic, Nonlinear Round-Key Generation for Block Ciphers"**
by Ayman E. Mohammed (manuscript under review).

A genetic algorithm (GA) is run **once, at design time**, to choose the parameters of GAKS, a deterministic,
invertible and nonlinear key-schedule family (crossover, two-way linear mixing, a partial S-box layer and
SHA-256-derived round constants). The resulting schedules are evaluated against the original DES and
PRESENT-80 key schedules.

## Requirements
Python 3.10+ with `numpy`, `scipy` (>= 1.11 for the MILP part) and `matplotlib`; a C compiler for `c_bench/`.

```
pip install -r requirements.txt
```

## Reproduce
```
python3 test_vectors.py            # DES and PRESENT-80 published test vectors
python3 nist.py                    # self-test of the SP 800-22 implementation against the spec examples
python3 ga_design.py des           # design-time GA (~5 min)  -> results/ga_des.json
python3 ga_design.py present       # (~6 min)                 -> results/ga_present.json
python3 evaluate.py                # experiments E1-E11, E14  -> results/eval.json
NSEQ=300 python3 run_nist.py DES-counter GAKS-DES-counter DES-random GAKS-DES-random PRESENT-counter GAKS-PRESENT-counter
python3 make_figures.py            # figures/ and the cost-matched ablation -> results/summary.json
python3 linear_bypass.py           # exact GF(2) S-box bypass length T_lin -> results/linear_bypass.json
python3 milp_rk.py 1,2,3,4,5,6,7 150 "PRESENT-80 single-key|PRESENT-80 related-key"
(cd c_bench && python3 gen_header.py && gcc -O2 -o bench bench.c && ./bench check && ./bench time)
```

The final GAKS-DES design is the GA chromosome in `results/ga_des.json` with the warm-up count fixed to
`w = 7` (`final_des_params()` in `gaks_core.py`); the decision was taken on a separate validation key set
(`results/validation_w_decision.json`).

**Seeds:** design keys 1, validation keys 777, evaluation keys 20261122, GA 2026.

## Official NIST STS 2.1.2
`run_nist.py` uses an independent Python implementation of nine SP 800-22 statistics. To run the
official reference implementation (all 15 tests), see the steps at the top of
`sts_official/run_official_sts.sh` (export streams, run `assess`, parse the reports).

## Files
| Path | Content |
| --- | --- |
| `gaks_core.py` | Vectorised DES, PRESENT-80, original and GAKS key schedules |
| `ga_design.py` | Design-time GA: chromosome, fitness function, operators |
| `evaluate.py` | Diffusion, SAC, related-key, complementation, weak-key, linearity and dependence experiments |
| `nist.py`, `run_nist.py` | SP 800-22 statistics (independent implementation) |
| `linear_bypass.py`, `milp_rk.py` | Exact S-box-bypass analysis and MILP models |
| `make_figures.py`, `gen_tables.py` | Figures, ablation and tables generated from the JSON results |
| `c_bench/` | C implementations of the four key schedules, cross-check and timing |
| `sts_official/` | Scripts for the official NIST STS 2.1.2 |
| `results/` | Raw results behind every table and figure |

## License and citation
MIT License (see `LICENSE`). If you use this code, please cite it using `CITATION.cff`.
