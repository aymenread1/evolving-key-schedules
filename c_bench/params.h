#include <stdint.h>
static const uint8_t PC1[56] = {57, 49, 41, 33, 25, 17, 9, 1, 58, 50, 42, 34, 26, 18, 10, 2, 59, 51, 43, 35, 27, 19, 11, 3, 60, 52, 44, 36, 63, 55, 47, 39, 31, 23, 15, 7, 62, 54, 46, 38, 30, 22, 14, 6, 61, 53, 45, 37, 29, 21, 13, 5, 28, 20, 12, 4};
static const uint8_t PC2[48] = {14, 17, 11, 24, 1, 5, 3, 28, 15, 6, 21, 10, 23, 19, 12, 4, 26, 8, 16, 7, 27, 20, 13, 2, 41, 52, 31, 37, 47, 55, 30, 40, 51, 45, 33, 48, 44, 49, 39, 56, 34, 53, 46, 42, 50, 36, 29, 32};
static const uint8_t SHIFTS[16] = {1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1};
static const uint8_t PSBOX[16] = {12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2};
#define DES_W 7
#define DES_STEPS 23
static const uint32_t DES_MASK[23] = {0xf84e901u, 0x7cd364du, 0x72c4749u, 0xc816e6eu, 0xcabd080u, 0x13bf8ffu, 0x8825dc5u, 0xc697e03u, 0x9ea779eu, 0x0a51791u, 0x323cf6du, 0xf71b382u, 0x6f55663u, 0x49e1861u, 0x80eaa5fu, 0x793ba6bu, 0xec0d48cu, 0x91ae1fbu, 0x0ce87d4u, 0x9eaea37u, 0x9733cb5u, 0x5cd4141u, 0x2832f89u};
static const uint8_t DES_ALPHA[23] = {3, 19, 20, 19, 23, 17, 5, 7, 25, 5, 9, 25, 5, 4, 17, 22, 3, 6, 14, 6, 9, 24, 4};
static const uint8_t DES_BETA[23] = {11, 11, 14, 1, 16, 3, 5, 3, 7, 15, 7, 4, 26, 22, 5, 2, 13, 11, 4, 6, 20, 5, 3};
#define DES_NNIB 2
static const uint8_t DES_NIB[2] = {11, 13};
static const uint64_t DES_RC[23] = {0x4da47cd3c1d842ull, 0x6773013a695b74ull, 0x2a73c2c1005e14ull, 0xffe787f8c06a77ull, 0x006418ecddc2e4ull, 0xf231f7ee4ff1b2ull, 0x2dbc4cd85d7a5dull, 0xfd56849630ae20ull, 0x3ceb73fb0614b6ull, 0xabad6511379176ull, 0x6cc5e2a73531c8ull, 0xc249c3048c6b90ull, 0x099ac696489408ull, 0xf936ea2a1737cbull, 0x5bacf022eee585ull, 0xbec8335fdef0adull, 0x55c68d8c44a474ull, 0x103510e615e712ull, 0xaf8606913eea49ull, 0x3a38f56f2071abull, 0x3796c4d9891e3bull, 0x77dce9b7027a17ull, 0xecb779a74dfb7aull};
#define P_W 7
#define P_NMASK 4
static const uint64_t P_MASK[4] = {0xbe3ced35a8ull, 0xfabb061e50ull, 0xd3fc7220bfull, 0x3ce448da4cull};
static const uint8_t P_ALPHA[4] = {4, 10, 2, 4};
static const uint8_t P_BETA[4] = {20, 9, 16, 20};
#define P_NNIBH 2
#define P_NNIBL 1
static const uint8_t P_NIBH[2] = {0, 9};
static const uint8_t P_NIBL[1] = {6};
static const uint64_t P_RCH[38] = {0x50df1e641aull, 0x99ae30e9a3ull, 0xd526939617ull, 0xfba8f36a1aull, 0x25101e4607ull, 0x504220484cull, 0xd8edea5389ull, 0x1ef80ad049ull, 0x8d9626f56dull, 0xcfe1cf2547ull, 0x4b8951aa8dull, 0x4de2bc3e09ull, 0x87fdfc2764ull, 0x835d247098ull, 0x3a8d494c81ull, 0xaee4b7f29eull, 0x1692f7122full, 0x69f1bd440full, 0x005aa944a7ull, 0x63cfe2e819ull, 0x3efe85120bull, 0xc37d273940ull, 0x0127f1076bull, 0x5fbbfe707cull, 0x2124d8a94dull, 0x55943d388aull, 0xeb69dee134ull, 0x70be4f64efull, 0x219028e7e5ull, 0xb8bdc3e4eeull, 0x24cd2f360full, 0xa45c774865ull, 0x9bcce322b1ull, 0xc7d4f3ffa7ull, 0xb88e8ead02ull, 0xd2282d369aull, 0x3c423722dbull, 0xa673da8d85ull};
static const uint64_t P_RCL[38] = {0x53336f5fa1ull, 0x5002428f8dull, 0xb0ae2f6fe6ull, 0x0437ac6ad6ull, 0xa912a97ad3ull, 0x610a6abc87ull, 0x4685c571b4ull, 0xbcc327515full, 0x938b2668c3ull, 0x6129203acfull, 0x596d3917d3ull, 0xcb85ccf2abull, 0xcce71e2e52ull, 0xe6c317bb33ull, 0x2700d969b6ull, 0x1bf4855b0bull, 0xfa1860c5c8ull, 0xb2535d8830ull, 0x784bef434eull, 0x005cfbfcf7ull, 0xf122d0c318ull, 0x5f26fc77c5ull, 0x2ad0b969e5ull, 0xaf7a08a60bull, 0xa646867457ull, 0x8477d48dd5ull, 0xc23cbd8cdbull, 0x1b338673a8ull, 0x127bfa3089ull, 0x49afff5590ull, 0xbcef93ec69ull, 0x009c1be6e4ull, 0x9fbc10d256ull, 0x997d15953bull, 0x050456c86full, 0x2f0302bcfbull, 0x550bec70ffull, 0x0baf820c1eull};
