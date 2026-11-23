"""
gift_core.py -- vectorised GIFT-64-128 (Banik et al., CHES 2017) with its original key schedule and the
GAKS-GIFT key schedule. Bit 0 = least significant bit, as in the GIFT specification.
The 128-bit key state K = k7||...||k0 (16-bit words) is held as two uint64 halves: H = k7..k4, L = k3..k0.
"""
import json
import numpy as np
from gaks_core import U, as_u64, popcount, _bits_to_int

GS = np.array([0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe], dtype=np.uint64)
GS_INV = np.argsort(GS).astype(np.uint64)
P64 = [0, 17, 34, 51, 48, 1, 18, 35, 32, 49, 2, 19, 16, 33, 50, 3,
       4, 21, 38, 55, 52, 5, 22, 39, 36, 53, 6, 23, 20, 37, 54, 7,
       8, 25, 42, 59, 56, 9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11,
       12, 29, 46, 63, 60, 13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15]
RCON = [0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3E, 0x3D, 0x3B, 0x37, 0x2F, 0x1E, 0x3C, 0x39, 0x33, 0x27, 0x0E,
        0x1D, 0x3A, 0x35, 0x2B, 0x16, 0x2C, 0x18, 0x30, 0x21, 0x02, 0x05, 0x0B]
ROUNDS = 28
M16 = U(0xFFFF); M32 = U(0xFFFFFFFF); M64 = U(0xFFFFFFFFFFFFFFFF)

def rotr16(x, s):
    s = U(s); return ((x >> s) | (x << (U(16) - s))) & M16

def rotl64(x, s):
    s = U(s % 64)
    return x if int(s) == 0 else ((x << s) | (x >> (U(64) - s))) & M64

def word(H, L, i):
    return ((L >> U(16 * i)) if i < 4 else (H >> U(16 * (i - 4)))) & M16

def gift_ks_update(H, L):
    """k7||...||k0 <- k1>>>2 || k0>>>12 || k7||...||k2"""
    k0, k1 = word(H, L, 0), word(H, L, 1)
    nL = (L >> U(32)) | ((H & M32) << U(32))          # k5..k2 move down to k3..k0
    nH = (H >> U(32)) | (rotr16(k0, 12) << U(32)) | (rotr16(k1, 2) << U(48))
    return nH, nL

def round_key(H, L):
    return L & M32                                     # U = k1 (bits 31..16), V = k0 (bits 15..0)

def ks_gift(kH, kL, nr=ROUNDS):
    H, L = as_u64(kH).copy(), as_u64(kL).copy()
    out = np.zeros((len(H), nr), dtype=np.uint64)
    for r in range(nr):
        out[:, r] = round_key(H, L)
        H, L = gift_ks_update(H, L)
    return out

# ------------------------------------------------------------------ cipher
_SPREAD = [(U(i), U(4 * i)) for i in range(16)]
def add_round_key(s, rk, r):
    u, v = rk >> U(16), rk & M16
    for i, p in _SPREAD:
        s = s ^ (((u >> i) & U(1)) << (p + U(1))) ^ (((v >> i) & U(1)) << p)
    c = RCON[r]
    s = s ^ (U(1) << U(63))
    for bit, pos in zip(range(6), [3, 7, 11, 15, 19, 23]):
        if (c >> bit) & 1:
            s = s ^ (U(1) << U(pos))
    return s

def subcells(s, box=GS):
    out = np.zeros_like(s)
    for j in range(16):
        sh = U(4 * j); out |= box[((s >> sh) & U(0xF)).astype(np.int64)] << sh
    return out

def permbits(s):
    out = np.zeros_like(s)
    for i, p in enumerate(P64):
        out |= ((s >> U(i)) & U(1)) << U(p)
    return out

def gift_encrypt(pt, rk, record=False, nr=ROUNDS):
    s = as_u64(pt).copy(); rec = []
    for r in range(nr):
        s = add_round_key(permbits(subcells(s)), rk[:, r], r)
        if record: rec.append(s.copy())
    return (s, np.stack(rec, 1)) if record else s

_PINV = [0] * 64
for i, p in enumerate(P64): _PINV[p] = i
def gift_decrypt(ct, rk, nr=ROUNDS):
    s = as_u64(ct).copy()
    for r in reversed(range(nr)):
        s = add_round_key(s, rk[:, r], r)               # XOR is an involution
        out = np.zeros_like(s)
        for p, i in enumerate(_PINV): out |= ((s >> U(p)) & U(1)) << U(i)
        s = subcells(out, GS_INV)
    return s

# ------------------------------------------------------------------ GAKS-GIFT
class GAKSGiftParams:
    """Chromosome: NBLOCKS gene blocks (64-bit crossover mask, two 6-bit rotation amounts 1..63) used
    cyclically, a 32-bit S-box mask over the 32 nibbles of the key state, and a 3-bit warm-up count.
    Round constants: the 5-bit step counter t+1 added to the low bits of both halves (hardware-friendly)."""
    MAXW = 7
    STEP = 64 + 6 + 6

    def __init__(self, masks, alpha, beta, smask, w, nblocks):
        self.masks, self.alpha, self.beta = list(masks), list(alpha), list(beta)
        self.smask, self.w, self.nblocks = int(smask), int(w), nblocks

    @staticmethod
    def nbits(nblocks):
        return nblocks * GAKSGiftParams.STEP + 32 + 3

    @classmethod
    def from_bits(cls, bits, nblocks=2):
        bits = [int(b) for b in bits]
        masks, alpha, beta = [], [], []
        for i in range(nblocks):
            g = bits[i * cls.STEP:(i + 1) * cls.STEP]
            masks.append(_bits_to_int(g[:64])); alpha.append(1 + _bits_to_int(g[64:70]) % 63); beta.append(1 + _bits_to_int(g[70:76]) % 63)
        off = nblocks * cls.STEP
        return cls(masks, alpha, beta, _bits_to_int(bits[off:off + 32]), _bits_to_int(bits[off + 32:off + 35]), nblocks)

    def sbox_count(self, nr=ROUNDS):
        return (self.w + nr - 1) * bin(self.smask).count("1")

    def to_dict(self):
        return {"nblocks": self.nblocks, "w": self.w, "smask": f"{self.smask:08x}", "masks": [f"{m:016x}" for m in self.masks],
                "alpha": self.alpha, "beta": self.beta, "sboxes_per_key_setup": self.sbox_count()}

def _sbox_half(x, sel):
    for j in range(16):
        if (sel >> j) & 1:
            sh = U(4 * j); x = (x & ~(U(0xF) << sh)) | (GS[((x >> sh) & U(0xF)).astype(np.int64)] << sh)
    return x

def gaks_gift_step(H, L, t, prm):
    H, L = gift_ks_update(H, L)                         # keep GIFT's word rotation
    g = t % prm.nblocks
    d = (H ^ L) & U(prm.masks[g]); H, L = H ^ d, L ^ d  # uniform crossover between halves
    L = L ^ rotl64(H, prm.alpha[g]); H = H ^ rotl64(L, prm.beta[g])   # two-way linear mixing
    H = _sbox_half(H, prm.smask >> 16); L = _sbox_half(L, prm.smask & 0xFFFF)  # partial S-box layer
    c = U(t + 1); H = H ^ c; L = L ^ c                   # counter constant
    return H, L

def ks_gaks_gift(kH, kL, prm, nr=ROUNDS):
    H, L = as_u64(kH).copy(), as_u64(kL).copy()
    out = np.zeros((len(H), nr), dtype=np.uint64); k = 0
    if prm.w == 0:
        out[:, 0] = round_key(H, L); k = 1
    t = 0
    while k < nr:
        H, L = gaks_gift_step(H, L, t, prm); t += 1
        if t >= prm.w:
            out[:, k] = round_key(H, L); k += 1
    return out

def key_halves(key128):
    return as_u64(key128 >> 64), as_u64(key128 & ((1 << 64) - 1))

if __name__ == "__main__":
    tv = [(0, 0, 0xf62bc3ef34f775ac),
          (0xfedcba9876543210fedcba9876543210, 0xfedcba9876543210, 0xc1b71f66160ff587)]
    for k, p, c in tv:
        H, L = key_halves(k)
        rk = ks_gift(H, L)
        ct = gift_encrypt(as_u64(p), rk)
        assert int(ct[0]) == c, hex(int(ct[0]))
        assert int(gift_decrypt(ct, rk)[0]) == p
    print("GIFT-64 official test vectors OK (giftcipher/gift, vectors 1 and 2)")
