#!/usr/bin/env bash
# Runs the official NIST STS 2.1.2 (all 15 tests, default parameters) on the six exported streams.
#
# 1. Download the reference code from NIST:
#      https://csrc.nist.gov/projects/random-bit-generation/documentation-and-software
#    (file sts-2_1_2.zip), unzip it, and build it:  cd sts-2.1.2/sts-2.1.2 && make
#    Create the output folders once (STS does not create them):
#      for t in Frequency BlockFrequency Runs LongestRun Rank FFT NonOverlappingTemplate \
#               OverlappingTemplate Universal LinearComplexity Serial ApproximateEntropy \
#               CumulativeSums RandomExcursions RandomExcursionsVariant; do
#        mkdir -p experiments/AlgorithmTesting/$t; done
# 2. From the repository root:  python3 sts_official/export_streams.py
# 3. Run:  bash sts_official/run_official_sts.sh /path/to/sts-2.1.2/sts-2.1.2
# 4. Parse: python3 sts_official/parse_sts.py   -> results/nist_official_sts.json
set -euo pipefail
STS_DIR="${1:?usage: run_official_sts.sh <path to the folder that contains ./assess>}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$HERE/reports"
for cfg in DES-counter GAKS-DES-counter DES-random GAKS-DES-random PRESENT-counter GAKS-PRESENT-counter; do
  echo "== $cfg"
  # answers to assess prompts: input file (0), path, all tests (1), no parameter changes (0),
  # 300 bitstreams, binary input format (1)
  ( cd "$STS_DIR" && printf "0\n%s\n1\n0\n300\n1\n" "$HERE/streams/$cfg.bin" | ./assess 1000000 > /dev/null )
  cp "$STS_DIR/experiments/AlgorithmTesting/finalAnalysisReport.txt" "$HERE/reports/$cfg.txt"
done
echo "Reports in $HERE/reports"
