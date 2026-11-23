"""gift_bench.py -- generates gift_bench.c (GIFT-64 original and GAKS-GIFT key schedules in C) from the Python
parameters, cross-checks the C round keys against Python, and times key set-up (10 runs). Writes gift_timing.json."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); os.chdir(ROOT)
import numpy as np
from gift_core import *

ev = json.load(open("results/gift/eval.json")); g = json.load(open("results/gift/ga_gift_nb2_s7.json"))
P = GAKSGiftParams.from_bits(g["chromosome"], g["nblocks"]); P.w = ev["final_w"]
nibH = [j for j in range(16) if (P.smask >> (16 + j)) & 1]; nibL = [j for j in range(16) if (P.smask >> j) & 1]
src = f"""#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
static const uint8_t GSB[16] = {{{', '.join(str(int(x)) for x in GS)}}};
static const uint64_t MASK[{P.nblocks}] = {{{', '.join(f'0x{m:016x}ull' for m in P.masks)}}};
static const int AL[{P.nblocks}] = {{{', '.join(map(str, P.alpha))}}}, BE[{P.nblocks}] = {{{', '.join(map(str, P.beta))}}};
static const int NIBH[{len(nibH)}] = {{{', '.join(map(str, nibH))}}}, NIBL[{len(nibL)}] = {{{', '.join(map(str, nibL))}}};
#define W {P.w}
#define NB {P.nblocks}
static inline uint16_t r16(uint16_t x, int s) {{ return (uint16_t)((x >> s) | (x << (16 - s))); }}
static inline uint64_t rl64(uint64_t x, int s) {{ return (x << s) | (x >> (64 - s)); }}
static inline void upd(uint64_t *H, uint64_t *L) {{
  uint16_t k0 = (uint16_t)(*L), k1 = (uint16_t)(*L >> 16);
  uint64_t nL = (*L >> 32) | (*H << 32);
  uint64_t nH = (*H >> 32) | ((uint64_t)r16(k0, 12) << 32) | ((uint64_t)r16(k1, 2) << 48);
  *H = nH; *L = nL;
}}
static void ks_gift(uint64_t H, uint64_t L, uint32_t rk[28]) {{
  for (int r = 0; r < 28; r++) {{ rk[r] = (uint32_t)L; upd(&H, &L); }}
}}
static inline uint64_t sbl(uint64_t x, const int *nib, int n) {{
  for (int i = 0; i < n; i++) {{ int s = 4 * nib[i]; x = (x & ~(0xFull << s)) | ((uint64_t)GSB[(x >> s) & 15] << s); }}
  return x;
}}
static void ks_gaks(uint64_t H, uint64_t L, uint32_t rk[28]) {{
  int k = 0, t = 0;
  if (W == 0) rk[k++] = (uint32_t)L;
  while (k < 28) {{
    upd(&H, &L); int g = t % NB;
    uint64_t d = (H ^ L) & MASK[g]; H ^= d; L ^= d;
    L ^= rl64(H, AL[g]); H ^= rl64(L, BE[g]);
    H = sbl(H, NIBH, {len(nibH)}); L = sbl(L, NIBL, {len(nibL)});
    H ^= (uint64_t)(t + 1); L ^= (uint64_t)(t + 1); t++;
    if (t >= W) rk[k++] = (uint32_t)L;
  }}
}}
static double now(void) {{ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }}
int main(int argc, char **argv) {{
  uint32_t rk[28];
  if (argc > 1 && !strcmp(argv[1], "check")) {{
    uint64_t keys[3][2] = {{{{0, 0}}, {{0xfedcba9876543210ull, 0xfedcba9876543210ull}}, {{0x0123456789abcdefull, 0x0f1e2d3c4b5a6978ull}}}};
    for (int a = 0; a < 2; a++) for (int i = 0; i < 3; i++) {{
      if (a) ks_gaks(keys[i][0], keys[i][1], rk); else ks_gift(keys[i][0], keys[i][1], rk);
      for (int r = 0; r < 28; r++) printf("%u%c", rk[r], r == 27 ? '\\n' : ',');
    }}
    return 0;
  }}
  const long N = 4000000; uint32_t acc = 0; const char *nm[2] = {{"GIFT-64", "GAKS-GIFT"}};
  printf("{{");
  for (int a = 0; a < 2; a++) {{
    double bn = 1e30, bc = 1e30;
    for (int rep = 0; rep < 5; rep++) {{
      uint64_t x = 0x9E3779B97F4A7C15ull; double t0 = now(); uint64_t c0 = __rdtsc();
      for (long i = 0; i < N; i++) {{
        x = x * 6364136223846793005ull + 1442695040888963407ull;
        if (a) ks_gaks(x, ~x, rk); else ks_gift(x, ~x, rk); acc ^= rk[27];
      }}
      double ns = (now() - t0) * 1e9 / N, cy = (double)(__rdtsc() - c0) / N;
      if (ns < bn) bn = ns;
      if (cy < bc) bc = cy;
    }}
    printf("%s\\"%s\\":{{\\"ns_per_key\\":%.1f,\\"tsc_cycles_per_key\\":%.0f}}", a ? "," : "", nm[a], bn, bc);
  }}
  printf(",\\"acc\\":%u}}\\n", acc);
  return 0;
}}
"""
open("c_bench/gift_bench.c", "w").write(src)
subprocess.run(["gcc", "-O2", "-w", "-o", "/tmp/gift_bench", "c_bench/gift_bench.c"], check=True)
out = subprocess.run(["/tmp/gift_bench", "check"], capture_output=True, text=True).stdout.split()
keys = [(0, 0), (0xfedcba9876543210, 0xfedcba9876543210), (0x0123456789abcdef, 0x0f1e2d3c4b5a6978)]
ref = [list(map(int, ks_gift(as_u64(h), as_u64(l))[0])) for h, l in keys] + [list(map(int, ks_gaks_gift(as_u64(h), as_u64(l), P)[0])) for h, l in keys]
got = [list(map(int, line.split(","))) for line in out]
print("C matches Python:", got == ref)
runs = [json.loads(subprocess.run(["/tmp/gift_bench"], capture_output=True, text=True).stdout) for _ in range(10)]
import statistics as st
res = {k: {"cycles_median": st.median(r[k]["tsc_cycles_per_key"] for r in runs), "ns_median": st.median(r[k]["ns_per_key"] for r in runs)} for k in ("GIFT-64", "GAKS-GIFT")}
ratios = [r["GAKS-GIFT"]["tsc_cycles_per_key"] / r["GIFT-64"]["tsc_cycles_per_key"] for r in runs]
res["ratio_range"] = [min(ratios), max(ratios)]; res["match"] = got == ref
json.dump(res, open("c_bench/gift_timing.json", "w"), indent=1); print(res)
