# Evolving Key Schedules Offline

Code, parameters, seeds and raw results for the article
**"Evolving Key Schedules Offline: A Cost-Aware Genetic Design of Deterministic, Nonlinear Round-Key Generation for Block Ciphers"**
by Ayman E. Mohammed (manuscript under review).

A genetic algorithm (GA) is run **once, at design time**, to choose the parameters of GAKS, a deterministic,
invertible and nonlinear key-schedule family (crossover, two-way linear mixing, a partial S-box layer and
step constants). The primary case study is **GIFT-64-128**, whose original key schedule is a bit permutation;
PRESENT-80 is a second case study, and DES is kept as a baseline (appendix of the article).

## Requirements
Python 3.10+ with `numpy`, `scipy` (>= 1.11 for the MILP part) and `matplotlib`; a C compiler for `c_bench/`.

```
pip install -r requirements.txt
```

## Reproduce: GIFT-64 (primary case study)
```
python3 gift_core.py               # GIFT-64-128 test vectors (designers' reference implementation)
SEED=7 POP=60 GENS=70 SUFFIX=_nb2_s7 python3 gift_design.py   # design-time GA -> results/gift/ga_gift_nb2_s7.json
python3 gift_evaluate.py nb2_s7    # validation rule (final w = 9) and E1-E11 -> results/gift/eval.json
python3 gift_linear.py             # exact GF(2) S-box bypass length -> results/gift/linear_bypass.json
python3 gift_ablation.py           # cost-matched ablation -> results/gift/ablation.json
pip install python-sat             # SAT model of related-key characteristics (CaDiCaL 1.9.5)
python3 gift_sat.py sk-orig 1,2,3,4,5,6,7,8,9           # single-key optima (designers' bounds)
python3 gift_sat.py rk-orig 1,2,3,4,5,6,7,8,9,10,11,12  # related-key, original schedule
python3 gift_sat.py rk-gaks-free 1,2,3,4,5,6,7,8,9 600  # GAKS-GIFT blocks starting where K1 is read
OFFSET=1 python3 gift_sat.py rk-gaks-free 8,9 600 6     # blocks starting one update later
python3 gift_sat.py ks-gaks 9 600 1                     # warm-up alone (lower bound)
python3 gift_sat_verify.py rk-gaks-free 8 11            # re-check a returned characteristic with gift_core.py
python3 gen_rk_table.py
python3 gift_milp.py 1,2,3,4,5,6 80; python3 gift_milp.py 7 450   # MILP cross-check (original schedule)
NSEQ=300 python3 run_nist.py GIFT-counter GAKSGIFT-counter GIFT-random GAKSGIFT-random
python3 c_bench/gift_bench.py      # C cross-check and timing -> c_bench/gift_timing.json
NANGATE_LIB=<path to NangateOpenCellLibrary_typical.lib> python3 hw/gift_hw.py   # Verilog, Yosys check, area
python3 make_gift_figures.py
```
Final GAKS-GIFT parameters: two gene blocks used alternately, masks `74324e75105437e8`, `97a2da2a17ba8f04`,
alpha = (61, 46), beta = (26, 11), S-box mask `0x00470080` (5 S-boxes per step), w = 9, step counter as constant;
180 S-box evaluations per key set-up.

Related-key bounds (minimum active S-boxes, `results/gift/sat_rk.json`): 7/8/9 rounds need 3/4/6 with the
original schedule and 6/11/>=13 with GAKS-GIFT (free-start blocks, warm-up excluded); full-cipher lower bound
43 versus 26 active S-boxes.

## Reproduce: DES and PRESENT-80
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

**Seeds:** design keys 1, validation keys 777, evaluation keys 20261122, GA 2026 (DES, PRESENT) and 7 (GIFT).

## Official NIST STS 2.1.2
The paper reports the official NIST STS 2.1.2 (all 15 tests, 188 sub-tests, 300 sequences of
10^6 bits per configuration). The raw reports are in `sts_official/reports/` and the summary in
`results/nist_official_sts.json`: GAKS-GIFT, GAKS-DES and GAKS-PRESENT pass 15/15 in every configuration;
GIFT-64 passes 0/15 (counter keys) and 1/15 (random keys), DES 1/15 and 5/15, PRESENT-80 1/15 (counter keys). To reproduce, follow the steps
at the top of `sts_official/run_official_sts.sh`. `run_nist.py` is an independent Python implementation
of nine of the statistics, used during design; it reproduces the official p-values for the
SP 800-22 reference input (binary expansion of e) to six decimals.

## Files
| Path | Content |
| --- | --- |
| `gift_core.py`, `gift_*.py`, `gen_rk_table.py` | GIFT-64-128, GAKS-GIFT, GA, evaluation, T_lin, MILP/SAT related-key bounds, ablation |
| `gaks_core.py` | Vectorised DES, PRESENT-80, original and GAKS key schedules |
| `ga_design.py` | Design-time GA: chromosome, fitness function, operators |
| `evaluate.py` | Diffusion, SAC, related-key, complementation, weak-key, linearity and dependence experiments |
| `nist.py`, `run_nist.py` | SP 800-22 statistics (independent implementation) |
| `linear_bypass.py`, `milp_rk.py` | Exact S-box-bypass analysis and MILP models |
| `make_figures.py`, `gen_tables.py` | Figures, ablation and tables generated from the JSON results |
| `c_bench/` | C implementations of the four key schedules, cross-check and timing (`timing_10runs.json`) |
| `sensitivity.py` | Fitness-weight sensitivity (GA reruns in `results/sensitivity/`) |
| `hw/` | Verilog generator, Yosys cross-check against Python, area with the NanGate 45 nm library (`synth_nangate.py`, needs `pip install pyosys` and the library file) |
| `evaluate_variant.py` | Hardware-friendly variant GAKS-DES-HW (two cyclic gene blocks, counter constants; `results/hw_variant/`) |
| `sts_official/` | Scripts for the official NIST STS 2.1.2 |
| `results/` | Raw results behind every table and figure |

## License and citation
MIT License (see `LICENSE`). If you use this code, please cite it using `CITATION.cff`.
