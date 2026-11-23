/* bench.c -- key-set-up cost of DES / GAKS-DES / PRESENT-80 / GAKS-PRESENT in portable C.
 * Bit permutations PC-1 and PC-2 use the same byte-wise lookup tables in both DES variants,
 * so the measured difference is the cost of the GAKS state update itself.
 * Build: gcc -O2 -o bench bench.c ; run: ./bench check | ./bench time */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <x86intrin.h>
#include "params.h"

static uint64_t pc1_t[8][256], pc2_t[7][256];

static uint64_t perm_slow(uint64_t x, const uint8_t *tab, int nout, int nin) {
    uint64_t o = 0;
    for (int k = 0; k < nout; k++) o |= ((x >> (nin - tab[k])) & 1ull) << (nout - 1 - k);
    return o;
}
static void init_tables(void) {
    for (int b = 0; b < 8; b++) for (int v = 0; v < 256; v++) pc1_t[b][v] = perm_slow((uint64_t)v << (8 * b), PC1, 56, 64);
    for (int b = 0; b < 7; b++) for (int v = 0; v < 256; v++) pc2_t[b][v] = perm_slow((uint64_t)v << (8 * b), PC2, 48, 56);
}
static inline uint64_t pc1(uint64_t k) { uint64_t o = 0; for (int b = 0; b < 8; b++) o |= pc1_t[b][(k >> (8 * b)) & 255]; return o; }
static inline uint64_t pc2(uint64_t s) { uint64_t o = 0; for (int b = 0; b < 7; b++) o |= pc2_t[b][(s >> (8 * b)) & 255]; return o; }
#define M28 0xFFFFFFFull
static inline uint64_t rotl28(uint64_t x, int s) { return ((x << s) | (x >> (28 - s))) & M28; }

static void ks_des(uint64_t key, uint64_t out[16]) {
    uint64_t s = pc1(key), C = s >> 28, D = s & M28;
    for (int r = 0; r < 16; r++) { C = rotl28(C, SHIFTS[r]); D = rotl28(D, SHIFTS[r]); out[r] = pc2((C << 28) | D); }
}
static void ks_gaks(uint64_t key, uint64_t out[16]) {
    uint64_t s = pc1(key), C = s >> 28, D = s & M28;
    for (int t = 0; t < DES_STEPS; t++) {
        int r = t - DES_W, sh = r >= 0 ? SHIFTS[r] : 1;
        C = rotl28(C, sh); D = rotl28(D, sh);
        uint64_t d = (C ^ D) & DES_MASK[t]; C ^= d; D ^= d;
        D ^= rotl28(C, DES_ALPHA[t]); C ^= rotl28(D, DES_BETA[t]);
        uint64_t S = (C << 28) | D;
        for (int i = 0; i < DES_NNIB; i++) { int sh4 = 4 * DES_NIB[i]; S = (S & ~(0xFull << sh4)) | ((uint64_t)PSBOX[(S >> sh4) & 15] << sh4); }
        S ^= DES_RC[t]; C = S >> 28; D = S & M28;
        if (r >= 0) out[r] = pc2(S);
    }
}
#define M40 0xFFFFFFFFFFull
static inline void rot80(uint64_t *H, uint64_t *L, int n) {
    n %= 80; if (n >= 40) { uint64_t t = *H; *H = *L; *L = t; n -= 40; }
    if (!n) return;
    uint64_t h = ((*H << n) | (*L >> (40 - n))) & M40, l = ((*L << n) | (*H >> (40 - n))) & M40; *H = h; *L = l;
}
static inline uint64_t rotl40(uint64_t x, int s) { return ((x << s) | (x >> (40 - s))) & M40; }
static void ks_present(uint64_t kh, uint64_t kl, uint64_t out[32]) {
    uint64_t H = kh, L = kl;
    for (int i = 1; i <= 32; i++) {
        out[i - 1] = (H << 24) | (L >> 16);
        if (i == 32) break;
        rot80(&H, &L, 61);
        H = (H & ~(0xFull << 36)) | ((uint64_t)PSBOX[(H >> 36) & 15] << 36);
        L ^= (uint64_t)i << 15;
    }
}
static inline uint64_t sbl(uint64_t x, const uint8_t *nib, int n) {
    for (int i = 0; i < n; i++) { int s = 4 * nib[i]; x = (x & ~(0xFull << s)) | ((uint64_t)PSBOX[(x >> s) & 15] << s); }
    return x;
}
static void ks_gpresent(uint64_t kh, uint64_t kl, uint64_t out[32]) {
    uint64_t H = kh, L = kl; int k = 0, total = P_W + 31;
    if (P_W == 0) out[k++] = (H << 24) | (L >> 16);
    for (int t = 0; t < total; t++) {
        rot80(&H, &L, 61);
        int g = t % P_NMASK; uint64_t d = (H ^ L) & P_MASK[g]; H ^= d; L ^= d;
        L ^= rotl40(H, P_ALPHA[g]); H ^= rotl40(L, P_BETA[g]);
        H = sbl(H, P_NIBH, P_NNIBH); L = sbl(L, P_NIBL, P_NNIBL);
        H ^= P_RCH[t]; L ^= P_RCL[t];
        if (t + 1 >= P_W && k < 32) out[k++] = (H << 24) | (L >> 16);
    }
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

int main(int argc, char **argv) {
    init_tables();
    if (argc > 1 && !strcmp(argv[1], "check")) {
        uint64_t keys[3] = {0x133457799BBCDFF1ull, 0x0123456789ABCDEFull, 0xFFFFFFFFFFFFFFFFull}, o[32];
        uint64_t pk[3][2] = {{0, 0}, {M40, M40}, {0x0123456789ull, 0xABCDEF0123ull}};
        printf("{\"des\":["); for (int i = 0; i < 3; i++) { ks_des(keys[i], o); printf("%s[", i ? "," : ""); for (int r = 0; r < 16; r++) printf("%s%llu", r ? "," : "", (unsigned long long)o[r]); printf("]"); }
        printf("],\"gaks\":["); for (int i = 0; i < 3; i++) { ks_gaks(keys[i], o); printf("%s[", i ? "," : ""); for (int r = 0; r < 16; r++) printf("%s%llu", r ? "," : "", (unsigned long long)o[r]); printf("]"); }
        printf("],\"present\":["); for (int i = 0; i < 3; i++) { ks_present(pk[i][0], pk[i][1], o); printf("%s[", i ? "," : ""); for (int r = 0; r < 32; r++) printf("%s%llu", r ? "," : "", (unsigned long long)o[r]); printf("]"); }
        printf("],\"gpresent\":["); for (int i = 0; i < 3; i++) { ks_gpresent(pk[i][0], pk[i][1], o); printf("%s[", i ? "," : ""); for (int r = 0; r < 32; r++) printf("%s%llu", r ? "," : "", (unsigned long long)o[r]); printf("]"); }
        printf("]}\n"); return 0;
    }
    const long N = 4000000; uint64_t o[32], acc = 0;
    const char *names[4] = {"DES", "GAKS-DES", "PRESENT-80", "GAKS-PRESENT"};
    printf("{");
    for (int a = 0; a < 4; a++) {
        double best_ns = 1e30, best_cyc = 1e30;
        for (int rep = 0; rep < 5; rep++) {
            uint64_t k = 0x9E3779B97F4A7C15ull; double t0 = now(); uint64_t c0 = __rdtsc();
            for (long i = 0; i < N; i++) {
                k = k * 6364136223846793005ull + 1442695040888963407ull;
                switch (a) { case 0: ks_des(k, o); break; case 1: ks_gaks(k, o); break;
                             case 2: ks_present(k >> 24, k & M40, o); break; default: ks_gpresent(k >> 24, k & M40, o); }
                acc ^= o[a < 2 ? 15 : 31];
            }
            double ns = (now() - t0) * 1e9 / N, cyc = (double)(__rdtsc() - c0) / N;
            if (ns < best_ns) best_ns = ns; if (cyc < best_cyc) best_cyc = cyc;
        }
        printf("%s\"%s\":{\"ns_per_key\":%.1f,\"tsc_cycles_per_key\":%.0f}", a ? "," : "", names[a], best_ns, best_cyc);
    }
    printf(",\"acc\":%llu}\n", (unsigned long long)acc);
    return 0;
}
