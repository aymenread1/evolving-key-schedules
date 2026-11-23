"""
gaks_core.py -- Vectorised DES / PRESENT-80 implementations with the standard
key schedules and the proposed GA-designed nonlinear key schedule (GAKS).

All bit numbering follows FIPS 46-3 (DES) and Bogdanov et al. CHES 2007 (PRESENT):
bit 1 of a DES table = most significant bit of the input word.
Every function accepts numpy uint64 arrays so thousands of keys are processed at once.
"""
import hashlib
import numpy as np

U = np.uint64
M28 = U((1 << 28) - 1)
M48 = U((1 << 48) - 1)
M56 = U((1 << 56) - 1)
M64 = U(0xFFFFFFFFFFFFFFFF)

# ---------------------------------------------------------------- DES tables
IP = [58,50,42,34,26,18,10,2,60,52,44,36,28,20,12,4,62,54,46,38,30,22,14,6,64,56,48,40,32,24,16,8,
      57,49,41,33,25,17,9,1,59,51,43,35,27,19,11,3,61,53,45,37,29,21,13,5,63,55,47,39,31,23,15,7]
FP = [40,8,48,16,56,24,64,32,39,7,47,15,55,23,63,31,38,6,46,14,54,22,62,30,37,5,45,13,53,21,61,29,
      36,4,44,12,52,20,60,28,35,3,43,11,51,19,59,27,34,2,42,10,50,18,58,26,33,1,41,9,49,17,57,25]
E = [32,1,2,3,4,5,4,5,6,7,8,9,8,9,10,11,12,13,12,13,14,15,16,17,16,17,18,19,20,21,20,21,
     22,23,24,25,24,25,26,27,28,29,28,29,30,31,32,1]
P = [16,7,20,21,29,12,28,17,1,15,23,26,5,18,31,10,2,8,24,14,32,27,3,9,19,13,30,6,22,11,4,25]
PC1 = [57,49,41,33,25,17,9,1,58,50,42,34,26,18,10,2,59,51,43,35,27,19,11,3,60,52,44,36,
       63,55,47,39,31,23,15,7,62,54,46,38,30,22,14,6,61,53,45,37,29,21,13,5,28,20,12,4]
PC2 = [14,17,11,24,1,5,3,28,15,6,21,10,23,19,12,4,26,8,16,7,27,20,13,2,
       41,52,31,37,47,55,30,40,51,45,33,48,44,49,39,56,34,53,46,42,50,36,29,32]
SHIFTS = [1,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1]
SBOX = [
 [14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7,0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8,
  4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0,15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13],
 [15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10,3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5,
  0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15,13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9],
 [10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8,13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1,
  13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7,1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12],
 [7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15,13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9,
  10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4,3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14],
 [2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9,14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6,
  4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14,11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3],
 [12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11,10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8,
  9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6,4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13],
 [4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1,13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6,
  1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2,6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12],
 [13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7,1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2,
  7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8,2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11]]
# reorganise S-boxes so that the 6-bit input is used directly as an index
SB = []
for s in SBOX:
    t = np.zeros(64, dtype=np.uint64)
    for x in range(64):
        row = ((x >> 4) & 2) | (x & 1)
        col = (x >> 1) & 0xF
        t[x] = s[row * 16 + col]
    SB.append(t)

PRESENT_S = np.array([0xC,5,6,0xB,9,0,0xA,0xD,3,0xE,0xF,8,4,7,1,2], dtype=np.uint64)
PRESENT_S_INV = np.argsort(PRESENT_S).astype(np.uint64)

def as_u64(x):
    return np.atleast_1d(np.asarray(x, dtype=np.uint64))

def perm(x, table, nin):
    """Generic bit permutation/selection, FIPS numbering (1 = MSB of an nin-bit word)."""
    x = as_u64(x)
    out = np.zeros_like(x)
    nout = len(table)
    for k, p in enumerate(table):
        out |= ((x >> U(nin - p)) & U(1)) << U(nout - 1 - k)
    return out

def rotl28(x, s):
    s = U(s)
    return ((x << s) | (x >> (U(28) - s))) & M28

def popcount(x):
    x = as_u64(x)
    b = x.view(np.uint8).reshape(x.shape + (8,))
    return np.unpackbits(b, axis=-1).sum(axis=-1).astype(np.int64)

# Non-parity (effective) key-bit positions of a 64-bit DES key, as bit index from LSB
EFFECTIVE_BITS = [63 - (p - 1) for p in PC1]  # the 56 bits PC-1 keeps

# ------------------------------------------------------------ DES key schedules
def ks_des(keys):
    k56 = perm(keys, PC1, 64)
    C, D = k56 >> U(28), k56 & M28
    out = np.zeros((len(k56), 16), dtype=np.uint64)
    for r in range(16):
        C, D = rotl28(C, SHIFTS[r]), rotl28(D, SHIFTS[r])
        out[:, r] = perm((C << U(28)) | D, PC2, 56)
    return out

def round_constants(n_rounds, label, bits):
    rc = []
    for t in range(1, n_rounds + 1):
        h = hashlib.sha256(f"{label}/RC/{t}".encode()).digest()
        rc.append(int.from_bytes(h, "big") >> (256 - bits))
    return rc

def _bits_to_int(bits):
    v = 0
    for b in bits:
        v = (v << 1) | int(b)
    return v

class GAKSParams:
    """Chromosome of the GA-designed DES key schedule (one gene block per state-update step).
    masks : (MAXW+16) 28-bit uniform-crossover masks
    alpha, beta : (MAXW+16) rotation amounts (1..27) of the two-way linear mixing
    smask : 14-bit mask selecting which nibbles of the 56-bit state pass through the S-box
    w     : number of warm-up (blank) steps before K1 is extracted (0..7)
    """
    MAXW = 7
    STEP = 28 + 5 + 5
    NBITS = (MAXW + 16) * STEP + 14 + 3

    def __init__(self, masks, alpha, beta, smask, w, rc=None, nblocks=None):
        self.masks, self.alpha, self.beta = list(masks), list(alpha), list(beta)
        self.smask, self.w = int(smask), int(w)
        self.rc = rc                      # None -> SHA-256 constants RC_DES; else a list of per-step constants
        self.nblocks = nblocks            # None -> one gene block per step; k -> k blocks used cyclically

    @staticmethod
    def nbits(nblocks=None):
        return (nblocks or (GAKSParams.MAXW + 16)) * GAKSParams.STEP + 14 + 3

    @classmethod
    def from_bits(cls, bits, nblocks=None, rc=None):
        bits = [int(b) for b in bits]
        nb = nblocks or (cls.MAXW + 16)
        blocks = []
        for i in range(nb):
            g = bits[i * cls.STEP:(i + 1) * cls.STEP]
            blocks.append((_bits_to_int(g[:28]), 1 + _bits_to_int(g[28:33]) % 27, 1 + _bits_to_int(g[33:38]) % 27))
        steps = [blocks[t % nb] for t in range(cls.MAXW + 16)]
        off = nb * cls.STEP
        smask = _bits_to_int(bits[off:off + 14])
        w = _bits_to_int(bits[off + 14:off + 17])
        return cls([b[0] for b in steps], [b[1] for b in steps], [b[2] for b in steps], smask, w, rc=rc, nblocks=nblocks)

    def sbox_count(self):
        return (self.w + 16) * bin(self.smask).count("1")

    def to_dict(self):
        n = self.w + 16
        return {"w": self.w, "smask": f"{self.smask:04x}",
                "masks": [f"{m:07x}" for m in self.masks[:n]],
                "alpha": self.alpha[:n], "beta": self.beta[:n],
                "sboxes_per_key_setup": self.sbox_count()}

W_FINAL_DES = 7   # final warm-up count: GA output (w=6) + 1, decided on a validation key set

def final_des_params(path="results/ga_des.json"):
    """Final GAKS-DES parameters: GA-selected genes with the warm-up count fixed to W_FINAL_DES."""
    import json as _json
    g = _json.load(open(path))
    p = GAKSParams.from_bits(g["chromosome"])
    p.w = W_FINAL_DES
    return p

RC_DES = round_constants(GAKSParams.MAXW + 16, "GAKS-DES", 56)
# hardware-friendly alternative: the 5-bit step counter t+1 added to the low bits of both halves
RC_CTR_DES = [((t + 1) << 28) | (t + 1) for t in range(GAKSParams.MAXW + 16)]

def sbox_layer(state, smask, nnib, sbox=PRESENT_S):
    for j in range(nnib):
        if (smask >> j) & 1:
            sh = U(4 * j)
            nib = (state >> sh) & U(0xF)
            state = (state & ~(U(0xF) << sh)) | (sbox[nib.astype(np.int64)] << sh)
    return state

def ks_gaks(keys, prm, use_rc=True):
    """Proposed deterministic, nonlinear key schedule for DES (16 x 48-bit subkeys)."""
    k56 = perm(keys, PC1, 64)
    C, D = k56 >> U(28), k56 & M28
    out = np.zeros((len(k56), 16), dtype=np.uint64)
    nrounds = prm.w + 16
    for t in range(nrounds):
        r = t - prm.w                     # main-round index (negative during warm-up)
        s = SHIFTS[r] if r >= 0 else 1
        C, D = rotl28(C, s), rotl28(D, s)
        m = U(prm.masks[t])               # GA operator 1: uniform crossover between halves
        diff = (C ^ D) & m
        C, D = C ^ diff, D ^ diff
        D = D ^ rotl28(C, prm.alpha[t])   # two-way linear mixing (invertible)
        C = C ^ rotl28(D, prm.beta[t])
        S = (C << U(28)) | D              # GA operator 2: nonlinear "mutation" (S-box layer)
        S = sbox_layer(S, prm.smask, 14)
        if use_rc:
            S ^= U(RC_DES[t] if prm.rc is None else prm.rc[t])   # round constant: breaks self-similarity
        C, D = S >> U(28), S & M28
        if r >= 0:
            out[:, r] = perm(S, PC2, 56)
    return out

# ------------------------------------------------------------ DES encryption
def des_encrypt(pt, subkeys, record=False):
    """pt: (n,) uint64; subkeys: (n,16) uint64. Returns ct (and per-round outputs if record)."""
    x = perm(pt, IP, 64)
    L, R = x >> U(32), x & U(0xFFFFFFFF)
    rounds = []
    for r in range(16):
        e = perm(R, E, 32) ^ subkeys[:, r]
        f = np.zeros_like(R)
        for i in range(8):
            six = ((e >> U(42 - 6 * i)) & U(0x3F)).astype(np.int64)
            f |= SB[i][six] << U(28 - 4 * i)
        f = perm(f, P, 32)
        L, R = R, L ^ f
        if record:
            rounds.append((L << U(32)) | R)
    ct = perm((R << U(32)) | L, FP, 64)
    if record:
        return ct, np.stack(rounds, axis=1)
    return ct

def des_decrypt(ct, subkeys):
    return des_encrypt(ct, subkeys[:, ::-1])

# ------------------------------------------------------------ PRESENT-80
PLAYER = [(16 * i) % 63 if i != 63 else 63 for i in range(64)]  # bit i -> PLAYER[i]

def present_sbox_layer(s):
    out = np.zeros_like(s)
    for j in range(16):
        sh = U(4 * j)
        out |= PRESENT_S[((s >> sh) & U(0xF)).astype(np.int64)] << sh
    return out

def present_player(s):
    out = np.zeros_like(s)
    for i in range(64):
        out |= ((s >> U(i)) & U(1)) << U(PLAYER[i])
    return out

M40 = U((1 << 40) - 1)

def _rot80_left(H, L, n):
    """80-bit register as two 40-bit halves (H = bits 79..40). Rotate left by n."""
    n %= 80
    if n >= 40:
        H, L = L, H
        n -= 40
    if n == 0:
        return H, L
    n = U(n)
    nH = ((H << n) | (L >> (U(40) - n))) & M40
    nL = ((L << n) | (H >> (U(40) - n))) & M40
    return nH, nL

def ks_present(kH, kL, nk=32):
    """Standard PRESENT-80 key schedule. kH,kL = upper/lower 40 bits. Returns (n,32) round keys."""
    H, L = as_u64(kH).copy(), as_u64(kL).copy()
    out = np.zeros((len(H), nk), dtype=np.uint64)
    for i in range(1, nk + 1):
        out[:, i - 1] = (H << U(24)) | (L >> U(16))
        if i == nk:
            break
        H, L = _rot80_left(H, L, 61)
        top = (H >> U(36)) & U(0xF)
        H = (H & ~(U(0xF) << U(36))) | (PRESENT_S[top.astype(np.int64)] << U(36))
        L ^= U(i) << U(15)               # bits 19..15 of the 80-bit register
    return out

def rotl40(x, s):
    s = U(s)
    return ((x << s) | (x >> (U(40) - s))) & M40

class GAKSPresentParams:
    """Same GA-designed structure instantiated on the 80-bit PRESENT key register.
    Gene blocks are cycled every NMASK steps."""
    NMASK = 4
    STEP = 40 + 6 + 6
    MAXW = 7
    NBITS = NMASK * STEP + 20 + 3

    def __init__(self, masks, alpha, beta, smask, w):
        self.masks, self.alpha, self.beta = list(masks), list(alpha), list(beta)
        self.smask, self.w = int(smask), int(w)

    @classmethod
    def from_bits(cls, bits):
        bits = [int(b) for b in bits]
        masks, alpha, beta = [], [], []
        for i in range(cls.NMASK):
            g = bits[i * cls.STEP:(i + 1) * cls.STEP]
            masks.append(_bits_to_int(g[:40]))
            alpha.append(1 + _bits_to_int(g[40:46]) % 39)
            beta.append(1 + _bits_to_int(g[46:52]) % 39)
        off = cls.NMASK * cls.STEP
        smask = _bits_to_int(bits[off:off + 20])
        w = _bits_to_int(bits[off + 20:off + 23])
        return cls(masks, alpha, beta, smask, w)

    def sbox_count(self, nk=32):
        return (self.w + nk - 1) * bin(self.smask).count("1")

    def to_dict(self):
        return {"w": self.w, "smask": f"{self.smask:05x}", "masks": [f"{m:010x}" for m in self.masks],
                "alpha": self.alpha, "beta": self.beta, "sboxes_per_key_setup": self.sbox_count()}

RC_PRESENT = round_constants(GAKSPresentParams.MAXW + 31, "GAKS-PRESENT", 80)

def ks_gaks_present(kH, kL, prm, nk=32):
    H, L = as_u64(kH).copy(), as_u64(kL).copy()
    out = np.zeros((len(H), nk), dtype=np.uint64)
    total = prm.w + nk - 1
    t_out = 0
    if prm.w == 0:
        out[:, 0] = (H << U(24)) | (L >> U(16)); t_out = 1
    for t in range(total):
        # (key j is extracted after step t with t + 1 >= w)
        H, L = _rot80_left(H, L, 61)
        g = t % prm.NMASK
        m = U(prm.masks[g])
        d = (H ^ L) & m
        H, L = H ^ d, L ^ d
        L = L ^ rotl40(H, prm.alpha[g])
        H = H ^ rotl40(L, prm.beta[g])
        H = sbox_layer(H, prm.smask >> 10, 10)
        L = sbox_layer(L, prm.smask & 0x3FF, 10)
        rc = RC_PRESENT[t]
        H ^= U(rc >> 40); L ^= U(rc & ((1 << 40) - 1))
        if t + 1 >= prm.w and t_out < nk:
            out[:, t_out] = (H << U(24)) | (L >> U(16)); t_out += 1
    return out

def present_encrypt(pt, rk, record=False, nrounds=31):
    s = as_u64(pt).copy()
    rounds = []
    for i in range(nrounds):
        s ^= rk[:, i]
        s = present_sbox_layer(s)
        s = present_player(s)
        if record:
            rounds.append(s.copy())          # state after round i+1, BEFORE the next round-key addition
    s ^= rk[:, nrounds]
    if record:
        return s, np.stack(rounds, axis=1)
    return s
