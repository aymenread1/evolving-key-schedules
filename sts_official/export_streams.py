"""
export_streams.py -- writes the six round-key streams used in the paper as packed binary files
(300 sequences x 1,000,000 bits each, most significant bit first) for the official NIST STS 2.1.2.
Usage (from the repository root):  python3 sts_official/export_streams.py
Output: sts_official/streams/<config>.bin  (37.5 MB each)
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); os.chdir(ROOT)
os.environ.setdefault("NSEQ", "300")
import run_nist

CONFIGS = ["DES-counter", "GAKS-DES-counter", "DES-random", "GAKS-DES-random", "PRESENT-counter", "GAKS-PRESENT-counter"]

def main(configs):
    out_dir = os.path.join(HERE, "streams"); os.makedirs(out_dir, exist_ok=True)
    for cfg in configs:
        path = os.path.join(out_dir, f"{cfg}.bin")
        with open(path, "wb") as f:
            for s in range(run_nist.NSEQ):
                bits = run_nist.stream(cfg, s)               # same streams as run_nist.py
                assert len(bits) == run_nist.NBITS
                f.write(np.packbits(bits).tobytes())         # MSB first, as STS binary mode expects
        print(cfg, os.path.getsize(path), "bytes")

if __name__ == "__main__":
    main(sys.argv[1:] or CONFIGS)
