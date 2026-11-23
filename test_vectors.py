"""Correctness checks against published test vectors."""
import numpy as np
from gaks_core import *

def h(x): return int(x)

# DES: classic worked example (key 133457799BBCDFF1)
k = as_u64(0x133457799BBCDFF1); p = as_u64(0x0123456789ABCDEF)
sk = ks_des(k)
c = des_encrypt(p, sk)
assert h(sk[0, 0]) == 0x1B02EFFC7072, hex(h(sk[0, 0]))
assert h(c[0]) == 0x85E813540F0AB405, hex(h(c[0]))
assert h(des_decrypt(c, sk)[0]) == 0x0123456789ABCDEF
# DES: NBS validation vector (key 0101..01, pt 95F8A5E5DD31D900 -> 8000000000000000)
c = des_encrypt(as_u64(0x95F8A5E5DD31D900), ks_des(as_u64(0x0101010101010101)))
assert h(c[0]) == 0x8000000000000000, hex(h(c[0]))
print("DES test vectors OK")

# PRESENT-80 test vectors (Bogdanov et al., CHES 2007, Appendix)
F40 = (1 << 40) - 1
vec = [(0, 0, 0, 0x5579C1387B228445), (F40, F40, 0, 0xE72C46C0F5945049),
       (0, 0, 0xFFFFFFFFFFFFFFFF, 0xA112FFC72F68417B), (F40, F40, 0xFFFFFFFFFFFFFFFF, 0x3333DCD3213210D2)]
for kh, kl, pt, ct in vec:
    rk = ks_present(as_u64(kh), as_u64(kl))
    out = present_encrypt(as_u64(pt), rk)
    assert h(out[0]) == ct, (hex(h(out[0])), hex(ct))
print("PRESENT-80 test vectors OK")
