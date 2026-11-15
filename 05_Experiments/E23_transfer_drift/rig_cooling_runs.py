"""Real rig, two runs under different cooling settings (L0_0 = baseline, L0_50 = current), 10 s logging.
Transfer-referenced drift per phase. Reference for phase W only: loss ratio from the oscilloscope
branch currents (RMS^2), 1.150 in L0_0 and 1.162 in L0_50 (E02), i.e. a drift of 1.010.
Development data (already inspected in E02). Pole and basis length are set for 10 s sampling."""
import sys, os
import numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
import tdrift

base_dir = sys.argv[1]
def run(name):
    x = pd.read_excel(os.path.join(base_dir, name, f"MAG1_1_{name}.xlsx"))
    x.columns = ["time"] + [str(c).split(" - ")[-1].strip() for c in x.columns[1:]]
    return x
B, N = run("L0_0"), run("L0_50")
COOL = {"U": "Tw_A-In", "V": "Tw_A-In", "W": "Tw_C-In"}
for a_s, M, skip in ((20.0, 4, 6), (30.0, 4, 6), (20.0, 3, 6)):
    tdrift.LOCK.update(a=float(np.exp(-10.0 / a_s)), M=M, skip=skip, n_jack=5, ma=3)
    print(f"Laguerre time constant {a_s:.0f} s, M {M}, skip {skip} samples")
    for ph in "UVW":
        names = [k for k in "1234" if f"{ph}1T{k}" in B and f"{ph}2T{k}" in B and f"{ph}1T{k}" in N]
        bp = [(f"T{k}", B[f"{ph}1T{k}"].to_numpy(float), B[f"{ph}2T{k}"].to_numpy(float)) for k in names]
        np_ = [(f"T{k}", N[f"{ph}1T{k}"].to_numpy(float), N[f"{ph}2T{k}"].to_numpy(float)) for k in names]
        r = tdrift.estimate_drift(B[COOL[ph]].to_numpy(float), bp, N[COOL[ph]].to_numpy(float), np_)
        lo, hi = r["ci95"]
        pos = " ".join(f"{np.exp(p['log_d']):.3f}" for p in r["positions"].values())
        ref = "  reference 1.010" if ph == "W" else ""
        print(f"  {ph}: drift {r['drift']:.3f} [{lo:.3f}, {hi:.3f}]  positions {pos}{ref}")
