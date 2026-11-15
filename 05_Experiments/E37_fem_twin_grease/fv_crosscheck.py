"""Finite-volume model (src/twin3d.py, uniform bottom cooling) with the corrected grease area resistance,
against the finite-element step responses of this folder (gap 40 mm, sensor below the grease).
python fv_crosscheck.py"""
import os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
import twin3d

Rg = 0.0225 * 0.086 * 0.070
for h in (1000, 8000):
    S = twin3d.step_response(twin3d.build(h, "bottom", dict(r_grease=Rg)), T_end=900)
    c_fv = S[1:, 4:].mean(1) / S[1:, :4].mean(1)
    F = np.loadtxt(os.path.join(here, f"fem_g40_h{h}_rb1_src1.csv"), delimiter=",")[:, :8]
    c_fe = F[1:, 4:].mean(1) / F[1:, :4].mean(1)
    print(f"h {h}: coupling at 10, 30, 120, 900 s  FV " + " ".join(f"{c_fv[t-1]:.4f}" for t in (10, 30, 120, 900))
          + "  FE " + " ".join(f"{c_fe[t-1]:.4f}" for t in (10, 30, 120, 900))
          + f"  | settled difference {abs(c_fv[-1] - c_fe[-1]):.4f}")
