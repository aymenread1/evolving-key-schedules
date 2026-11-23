import os, sys, re, subprocess, numpy as np
sys.path.insert(0, "hw"); sys.path.insert(0, ".")
import gen_verilog as GV
from gaks_core import *
prm, nb = GV.load_hw("hw2")
rng = np.random.default_rng(11); M28i=(1<<28)-1
rot=lambda x,s,w: ((x<<s)|(x>>(w-s)))&((1<<w)-1); SBi=[int(v) for v in PRESENT_S]
def ref(S,t):
    C,D=S>>28,S&M28i; r=t-prm.w; s=SHIFTS[r] if r>=0 else 1
    C,D=rot(C,s,28),rot(D,s,28); d=(C^D)&prm.masks[t]; C^=d; D^=d
    D^=rot(C,prm.alpha[t],28); C^=rot(D,prm.beta[t],28); S=(C<<28)|D
    for j in range(14):
        if (prm.smask>>j)&1: S=(S&~(0xF<<4*j))|(SBi[(S>>4*j)&0xF]<<4*j)
    return S ^ prm.rc[t]
cases=[(int(rng.integers(0,2**56,dtype=np.uint64)),t) for t in range(prm.w+16)]
cmds=["read_verilog hw/ks_gaks_des_hw.v","synth -flatten -noabc -top gaks_des_hw_step"]+[f"eval -set si 56'h{S:014x} -set t 5'd{t} -show so" for S,t in cases]
out=subprocess.run(["yowasp-yosys","-p","; ".join(cmds)],capture_output=True,text=True).stdout
got=[int(m,2) for m in re.findall(r"\\so = \d+'([01]+)",out)]
print(len(got), all(a==ref(*c) for a,c in zip(got,cases)) and len(got)==len(cases))
# also ks_gaks consistency: python step == ks_gaks first round key
k=0x133457799BBCDFF1; S=int(perm(as_u64(k),PC1,64)[0])
for t in range(prm.w+1): S=ref(S,t)
print(int(perm(as_u64(S),PC2,56)[0])==int(ks_gaks(as_u64(k),prm)[0,0]))
