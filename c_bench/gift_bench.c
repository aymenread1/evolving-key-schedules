#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
static const uint8_t GSB[16] = {1, 10, 4, 12, 6, 15, 3, 9, 2, 13, 11, 7, 5, 0, 8, 14};
static const uint64_t MASK[2] = {0x74324e75105437e8ull, 0x97a2da2a17ba8f04ull};
static const int AL[2] = {61, 46}, BE[2] = {26, 11};
static const int NIBH[4] = {0, 1, 2, 6}, NIBL[1] = {7};
#define W 9
#define NB 2
static inline uint16_t r16(uint16_t x, int s) { return (uint16_t)((x >> s) | (x << (16 - s))); }
static inline uint64_t rl64(uint64_t x, int s) { return (x << s) | (x >> (64 - s)); }
static inline void upd(uint64_t *H, uint64_t *L) {
  uint16_t k0 = (uint16_t)(*L), k1 = (uint16_t)(*L >> 16);
  uint64_t nL = (*L >> 32) | (*H << 32);
  uint64_t nH = (*H >> 32) | ((uint64_t)r16(k0, 12) << 32) | ((uint64_t)r16(k1, 2) << 48);
  *H = nH; *L = nL;
}
static void ks_gift(uint64_t H, uint64_t L, uint32_t rk[28]) {
  for (int r = 0; r < 28; r++) { rk[r] = (uint32_t)L; upd(&H, &L); }
}
static inline uint64_t sbl(uint64_t x, const int *nib, int n) {
  for (int i = 0; i < n; i++) { int s = 4 * nib[i]; x = (x & ~(0xFull << s)) | ((uint64_t)GSB[(x >> s) & 15] << s); }
  return x;
}
static void ks_gaks(uint64_t H, uint64_t L, uint32_t rk[28]) {
  int k = 0, t = 0;
  if (W == 0) rk[k++] = (uint32_t)L;
  while (k < 28) {
    upd(&H, &L); int g = t % NB;
    uint64_t d = (H ^ L) & MASK[g]; H ^= d; L ^= d;
    L ^= rl64(H, AL[g]); H ^= rl64(L, BE[g]);
    H = sbl(H, NIBH, 4); L = sbl(L, NIBL, 1);
    H ^= (uint64_t)(t + 1); L ^= (uint64_t)(t + 1); t++;
    if (t >= W) rk[k++] = (uint32_t)L;
  }
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
int main(int argc, char **argv) {
  uint32_t rk[28];
  if (argc > 1 && !strcmp(argv[1], "check")) {
    uint64_t keys[3][2] = {{0, 0}, {0xfedcba9876543210ull, 0xfedcba9876543210ull}, {0x0123456789abcdefull, 0x0f1e2d3c4b5a6978ull}};
    for (int a = 0; a < 2; a++) for (int i = 0; i < 3; i++) {
      if (a) ks_gaks(keys[i][0], keys[i][1], rk); else ks_gift(keys[i][0], keys[i][1], rk);
      for (int r = 0; r < 28; r++) printf("%u%c", rk[r], r == 27 ? '\n' : ',');
    }
    return 0;
  }
  const long N = 4000000; uint32_t acc = 0; const char *nm[2] = {"GIFT-64", "GAKS-GIFT"};
  printf("{");
  for (int a = 0; a < 2; a++) {
    double bn = 1e30, bc = 1e30;
    for (int rep = 0; rep < 5; rep++) {
      uint64_t x = 0x9E3779B97F4A7C15ull; double t0 = now(); uint64_t c0 = __rdtsc();
      for (long i = 0; i < N; i++) {
        x = x * 6364136223846793005ull + 1442695040888963407ull;
        if (a) ks_gaks(x, ~x, rk); else ks_gift(x, ~x, rk); acc ^= rk[27];
      }
      double ns = (now() - t0) * 1e9 / N, cy = (double)(__rdtsc() - c0) / N;
      if (ns < bn) bn = ns;
      if (cy < bc) bc = cy;
    }
    printf("%s\"%s\":{\"ns_per_key\":%.1f,\"tsc_cycles_per_key\":%.0f}", a ? "," : "", nm[a], bn, bc);
  }
  printf(",\"acc\":%u}\n", acc);
  return 0;
}
