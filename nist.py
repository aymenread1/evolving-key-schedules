"""
nist.py -- Independent numpy implementation of nine NIST SP 800-22 Rev.1a tests
(validated against the worked examples of the specification in selftest()).
"""
import numpy as np
from scipy.special import erfc, gammaincc
from scipy.stats import norm

def frequency(e):
    n = len(e); s = np.abs(np.sum(2 * e.astype(np.int64) - 1)) / np.sqrt(n)
    return erfc(s / np.sqrt(2))

def block_frequency(e, M=128):
    n = len(e); N = n // M
    pi = e[:N * M].reshape(N, M).mean(axis=1)
    chi = 4 * M * np.sum((pi - 0.5) ** 2)
    return gammaincc(N / 2, chi / 2)

def cusum(e, mode=0):
    x = 2 * e.astype(np.int64) - 1
    if mode == 1: x = x[::-1]
    S = np.cumsum(x); z = np.max(np.abs(S)); n = len(e)
    k1 = np.arange(int((-n / z + 1) // 4), int((n / z - 1) // 4) + 1)
    k2 = np.arange(int((-n / z - 3) // 4), int((n / z - 1) // 4) + 1)
    sq = np.sqrt(n)
    p = 1 - np.sum(norm.cdf((4 * k1 + 1) * z / sq) - norm.cdf((4 * k1 - 1) * z / sq)) \
          + np.sum(norm.cdf((4 * k2 + 3) * z / sq) - norm.cdf((4 * k2 + 1) * z / sq))
    return float(p)

def runs(e):
    n = len(e); pi = e.mean()
    if abs(pi - 0.5) >= 2 / np.sqrt(n): return 0.0
    v = 1 + np.count_nonzero(e[1:] != e[:-1])
    return erfc(abs(v - 2 * n * pi * (1 - pi)) / (2 * np.sqrt(2 * n) * pi * (1 - pi)))

def longest_run(e):
    n = len(e)
    if n < 6272: M, lo, hi, pi = 8, 1, 4, [0.21484375, 0.3671875, 0.23046875, 0.1875]
    elif n < 750000: M, lo, hi, pi = 128, 4, 9, [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
    else: M, lo, hi, pi = 10000, 10, 16, [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    N = n // M; blocks = e[:N * M].reshape(N, M).astype(np.int8)
    # longest run of ones per block (vectorised)
    pad = np.zeros((N, 1), np.int8)
    d = np.diff(np.concatenate([pad, blocks, pad], axis=1), axis=1)
    longest = np.zeros(N, dtype=np.int64)
    for i in range(N):
        st = np.flatnonzero(d[i] == 1); en = np.flatnonzero(d[i] == -1)
        longest[i] = (en - st).max() if len(st) else 0
    v = np.bincount(np.clip(longest, lo, hi) - lo, minlength=hi - lo + 1)
    pi = np.array(pi); chi = np.sum((v - N * pi) ** 2 / (N * pi))
    return gammaincc((hi - lo) / 2, chi / 2)

def dft(e):
    n = len(e); x = 2 * e.astype(np.float64) - 1
    m = np.abs(np.fft.fft(x))[: n // 2]
    T = np.sqrt(np.log(1 / 0.05) * n); N0 = 0.95 * n / 2; N1 = np.count_nonzero(m < T)
    d = (N1 - N0) / np.sqrt(n * 0.95 * 0.05 / 4)
    return erfc(abs(d) / np.sqrt(2))

def _pattern_counts(e, m):
    if m == 0: return np.array([len(e)])
    n = len(e); ext = np.concatenate([e, e[:m - 1]]).astype(np.int64)
    idx = np.zeros(n, dtype=np.int64)
    for j in range(m):
        idx = (idx << 1) | ext[j:j + n]
    return np.bincount(idx, minlength=2 ** m)

def approximate_entropy(e, m=10):
    n = len(e)
    def phi(mm):
        c = _pattern_counts(e, mm) / n; c = c[c > 0]
        return np.sum(c * np.log(c))
    ap = phi(m) - phi(m + 1)
    chi = 2 * n * (np.log(2) - ap)
    return gammaincc(2 ** (m - 1), chi / 2)

def serial(e, m=16):
    n = len(e)
    def psi(mm):
        if mm <= 0: return 0.0
        c = _pattern_counts(e, mm)
        return (2 ** mm) / n * np.sum(c.astype(np.float64) ** 2) - n
    p0, p1, p2 = psi(m), psi(m - 1), psi(m - 2)
    d1 = p0 - p1; d2 = p0 - 2 * p1 + p2
    return gammaincc(2 ** (m - 2), d1 / 2), gammaincc(2 ** (m - 3), d2 / 2)

# NOTE: the DFT (spectral) test is implemented but excluded from the battery because the
# worked examples in SP 800-22 Rev.1a are not mutually consistent for it (see README).
TESTS = ["Frequency", "BlockFrequency", "CumulativeSums(fwd)", "CumulativeSums(rev)", "Runs",
         "LongestRun", "ApproximateEntropy", "Serial-1", "Serial-2"]

def run_all(e):
    s1, s2 = serial(e)
    return [frequency(e), block_frequency(e), cusum(e, 0), cusum(e, 1), runs(e), longest_run(e),
            approximate_entropy(e), s1, s2]

def selftest():
    b = lambda s: np.array([int(c) for c in s], dtype=np.uint8)
    eps = b("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000")
    chk = {"frequency": (frequency(eps), 0.109599), "block_frequency(M=10)": (block_frequency(eps, 10), 0.706438),
           "cusum fwd": (cusum(eps, 0), 0.219194), "cusum rev": (cusum(eps, 1), 0.114866),
           "runs": (runs(eps), 0.500798),
           "apen(m=2)": (approximate_entropy(eps, 2), 0.235301)}
    lr = b("11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010")
    chk["longest_run"] = (longest_run(lr), 0.180609)
    s1, s2 = serial(b("0011011101"), 3)
    chk["serial p1"] = (s1, 0.808792); chk["serial p2"] = (s2, 0.670320)
    ok = True
    for k, (got, exp) in chk.items():
        good = abs(got - exp) < 2e-3
        ok &= good
        print(f"{k:24s} got {got:.6f} expected {exp:.6f} {'OK' if good else 'MISMATCH'}")
    return ok

if __name__ == "__main__":
    print("all OK" if selftest() else "SOME MISMATCH")
