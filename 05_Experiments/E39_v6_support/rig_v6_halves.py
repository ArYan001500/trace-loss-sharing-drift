"""Frozen method v6 on the real mission recording: baseline = first half, current = second half, per phase.
No change of the loss ratio is expected within one 10-minute run, so the interval should contain 1.
Also reported: the same with the halves swapped, and the absolute two-temperature ratios of each half
(transfer drift needs no current; v6 guards use the phase current shifted by the clock offset of 11 s)."""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from load_mission import load
import method_v6

df = load(sys.argv[1], 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
CIRC = {"U": "A", "V": "A", "W": "C"}
n = len(df); h = n // 2
for ph in "UVW":
    I = df[PH[ph]].shift(-11).ffill().bfill().to_numpy()
    Tc = (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
    pairs = [(f"T{k}", df[f"{ph}1T{k}"].to_numpy(), df[f"{ph}2T{k}"].to_numpy()) for k in "1234"
             if f"{ph}1T{k}" in df and f"{ph}2T{k}" in df]
    cut = lambda s: [(nm, a[s], b[s]) for nm, a, b in pairs]
    for label, sb, sn in (("first->second", slice(0, h), slice(h, n)), ("second->first", slice(h, n), slice(0, h))):
        r = method_v6.estimate_drift((I[sb], Tc[sb], cut(sb)), (I[sn], Tc[sn], cut(sn)))
        if "drift" in r:
            lo, hi = r["ci95"]
            print(f"{ph} {label}: drift {r['drift']:.4f} [{lo:.4f}, {hi:.4f}] contains 1: {lo <= 1 <= hi}; v3 {r['v3']:.4f}; "
                  f"z_x {r['zx']:+.2f}, shape z {r['shape_z']:+.2f}; {'ABSTAIN: ' + r['reason'] if r['abstain'] else 'reported'}")
        else:
            print(f"{ph} {label}: ABSTAIN: {r.get('reason')}")
