"""Sensitivity of the rig result (first half -> second half drift) to analysis choices.
Transfer drift: Laguerre time constant 5/10/20 s, basis length 4/6/8, skipped samples 30/60/90 (27 settings).
Full v6 decision: electrical clock offset 10/11/12 s (affects only the guards, which use the current).
The frozen settings (10 s, 6, 60, offset 11) are the reference; nothing here changes the method."""
import sys, os, itertools
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from load_mission import load
import tdrift, method_v6

df = load(sys.argv[1], 0); n = len(df); h = n // 2
CIRC = {"U": "A", "V": "A", "W": "C"}; PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
frozen = dict(tdrift.LOCK)
for ph in "UVW":
    Tc = (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
    pairs = [(f"T{k}", df[f"{ph}1T{k}"].to_numpy(), df[f"{ph}2T{k}"].to_numpy()) for k in "1234" if f"{ph}1T{k}" in df and f"{ph}2T{k}" in df]
    cut = lambda s: [(nm, a[s], b[s]) for nm, a, b in pairs]
    ds = []
    for ts, M, sk in itertools.product((5, 10, 20), (4, 6, 8), (30, 60, 90)):
        tdrift.LOCK.update(a=float(np.exp(-1 / ts)), M=M, skip=sk)
        ds.append(tdrift.estimate_drift(Tc[:h], cut(slice(0, h)), Tc[h:], cut(slice(h, n)))["drift"])
    tdrift.LOCK.update(frozen)
    ds = np.array(ds)
    print(f"{ph}: transfer drift over 27 settings: min {ds.min():.4f}, median {np.median(ds):.4f}, max {ds.max():.4f}")
    for off in (10, 11, 12):
        I = df[PH[ph]].shift(-off).ffill().bfill().to_numpy()
        r = method_v6.estimate_drift((I[:h], Tc[:h], cut(slice(0, h))), (I[h:], Tc[h:], cut(slice(h, n))))
        print(f"   clock offset {off} s: drift {r['drift']:.4f} [{r['ci95'][0]:.4f}, {r['ci95'][1]:.4f}], "
              f"z_x {r['zx']:+.2f}, shape z {r['shape_z']:+.2f}, {'ABSTAIN' if r['abstain'] else 'reported'}")
