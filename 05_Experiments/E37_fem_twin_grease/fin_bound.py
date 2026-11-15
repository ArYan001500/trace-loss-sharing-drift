"""Closed-form bound for the SETTLED thermal coupling between neighbouring modules on a cooled plate.
A plate of conductivity k and thickness t, cooled on one face with heat-transfer coefficient h, conducts heat
sideways like a fin: the steady excess temperature decays as exp(-x / lambda), lambda = sqrt(k t / h).
The neighbour's footprint edge is at most as warm as the centre of its own footprint, which by symmetry equals
the own response at the sensor. Two-dimensional spreading only adds decay. Hence
    c_settled <= exp(-d / lambda),    d = distance from the sensor to the nearest edge of the neighbour footprint.
Checked against the settled (900 s) coupling of every finite-element variant (E37, corrected grease; copy of E33).
python fin_bound.py FEM_DIR"""
import sys, os, re
import numpy as np

k, t = 167.0, 0.025
fem = sys.argv[1]
print(f"{'file':34s} {'gap mm':>6s} {'h':>6s} {'lambda mm':>9s} {'FEM settled':>11s} {'fin bound':>9s} {'bound holds':>11s}")
for fn in sorted(os.listdir(fem)):
    m = re.match(r"fem_(?:g(\d+)_)?h(\d+)_rb1_src1\.csv$", fn) or re.match(r"fem_step_h(\d+)\.csv$", fn)
    if not m:
        continue
    if fn.startswith("fem_step"):
        gap, h = 40, int(m.group(1))
    else:
        gap, h = int(m.group(1) or 40), int(m.group(2))
    S = np.loadtxt(os.path.join(fem, fn), delimiter=",")[:, :8]
    c = S[-1, 4:].mean() / S[-1, :4].mean()
    d = 0.043 + gap / 1000; lam = np.sqrt(k * t / h); b = np.exp(-d / lam)
    print(f"{fn:34s} {gap:6d} {h:6d} {1e3*lam:9.1f} {c:11.4f} {b:9.4f} {str(c <= b):>11s}")
print("\nInterval term of TRACE with this bound: |log d| * 2 c / (1 - c). For the rig (gap 40 mm):")
for h in (500, 1000, 2000, 8000):
    lam = np.sqrt(k * t / h); b = np.exp(-0.083 / lam)
    print(f"   h {h:5d}: c_max {b:.3f} -> extension {2*b/(1-b):.3f} x |log d|")
