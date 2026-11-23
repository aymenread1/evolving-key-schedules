"""Generates params.h for the C benchmark from the Python reference (no hand transcription)."""
import sys, json
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
import os; os.chdir("..")
from gaks_core import *

P = final_des_params()
G = json.load(open("results/ga_present.json")); PP = GAKSPresentParams.from_bits(G["chromosome"])
L = []
arr = lambda name, typ, vals: L.append(f"static const {typ} {name}[{len(vals)}] = {{{', '.join(map(str, vals))}}};")
arr("PC1", "uint8_t", PC1); arr("PC2", "uint8_t", PC2); arr("SHIFTS", "uint8_t", SHIFTS)
arr("PSBOX", "uint8_t", [int(x) for x in PRESENT_S])
n = P.w + 16
L.append(f"#define DES_W {P.w}\n#define DES_STEPS {n}")
arr("DES_MASK", "uint32_t", [f"0x{m:07x}u" for m in P.masks[:n]])
arr("DES_ALPHA", "uint8_t", P.alpha[:n]); arr("DES_BETA", "uint8_t", P.beta[:n])
nib = [j for j in range(14) if (P.smask >> j) & 1]
L.append(f"#define DES_NNIB {len(nib)}"); arr("DES_NIB", "uint8_t", nib)
arr("DES_RC", "uint64_t", [f"0x{r:014x}ull" for r in RC_DES[:n]])
L.append(f"#define P_W {PP.w}\n#define P_NMASK {PP.NMASK}")
arr("P_MASK", "uint64_t", [f"0x{m:010x}ull" for m in PP.masks]); arr("P_ALPHA", "uint8_t", PP.alpha); arr("P_BETA", "uint8_t", PP.beta)
nibH = [j for j in range(10) if (PP.smask >> (10 + j)) & 1]; nibL = [j for j in range(10) if (PP.smask >> j) & 1]
L.append(f"#define P_NNIBH {len(nibH)}\n#define P_NNIBL {len(nibL)}")
arr("P_NIBH", "uint8_t", nibH or [0]); arr("P_NIBL", "uint8_t", nibL or [0])
rc = RC_PRESENT[:PP.w + 31]
arr("P_RCH", "uint64_t", [f"0x{r >> 40:010x}ull" for r in rc]); arr("P_RCL", "uint64_t", [f"0x{r & ((1 << 40) - 1):010x}ull" for r in rc])
open("c_bench/params.h", "w").write("#include <stdint.h>\n" + "\n".join(L) + "\n")

# reference outputs for cross-checking
keys = [0x133457799BBCDFF1, 0x0123456789ABCDEF, 0xFFFFFFFFFFFFFFFF]
ref = {"des": [[int(x) for x in ks_des(as_u64(k))[0]] for k in keys],
       "gaks": [[int(x) for x in ks_gaks(as_u64(k), P)[0]] for k in keys]}
pk = [(0, 0), ((1 << 40) - 1, (1 << 40) - 1), (0x0123456789, 0xABCDEF0123)]
ref["present"] = [[int(x) for x in ks_present(as_u64(h), as_u64(l))[0]] for h, l in pk]
ref["gpresent"] = [[int(x) for x in ks_gaks_present(as_u64(h), as_u64(l), PP)[0]] for h, l in pk]
json.dump(ref, open("c_bench/reference.json", "w"))
print("ok", P.w, nib, PP.w, nibH, nibL)
