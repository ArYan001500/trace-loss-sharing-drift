"""Closed-form screening estimate of thermal coupling between neighbouring heat sources, checked against the
finite-element twin (E37, corrected grease; copy of E33) and the real two-transistor board (E25).

Derivation (1D diffusion from a plane source of constant power switched on at t = 0, in a medium of
diffusivity alpha; Carslaw and Jaeger, conduction of heat in solids):
    T(x, t) = (2 q / k) sqrt(alpha t) ierfc(eta),   eta = x / (2 sqrt(alpha t)),
    ierfc(eta) = exp(-eta^2) / sqrt(pi) - eta erfc(eta).
Slope at distance d over slope at the source (dT/dt = (q/k) sqrt(alpha/(pi t)) exp(-eta^2)):
    c_slope(t) = exp(-eta^2) = exp(-1/(4 Fo)),   Fo = alpha t / d^2.
    (E33 used erfc(eta) here, which is the heat-flux ratio, not the slope ratio; corrected in this copy.)
Rise at distance d over rise at the source:
    c_rise(t) = sqrt(pi) ierfc(eta) = exp(-eta^2) - sqrt(pi) eta erfc(eta).
Three-dimensional spreading and heat removal by the coolant only lower the neighbour's response (comparison
principle for the heat equation with non-negative sources and sinks), while the own response at early times is
set by the local heat capacity. The 1D values are screening estimates, not bounds for real hardware.
Screening rule: c_slope(w) <= eps  <=>  Fo = alpha w / d^2 <= 1 / (4 ln(1/eps)).
python coupling_bound.py FEM_DIR SANZ_NPZ"""
import sys, os
import numpy as np
from scipy.special import erfc

def c_slope(t, d, a): return np.exp(-d * d / (4 * a * t))
def c_rise(t, d, a):
    e = d / (2 * np.sqrt(a * t)); return np.exp(-e * e) - np.sqrt(np.pi) * e * erfc(e)

a_al = 167 / (2700 * 900)
print(f"aluminium diffusivity {a_al:.3e} m2/s")
for eps in (0.01, 0.05):
    print(f"screening rule for slope coupling <= {eps}: Fo = alpha w / d^2 <= {1 / (4 * np.log(1 / eps)):.4f}")

fem = sys.argv[1]
ts = np.array([4, 10, 20, 30, 60])
print("\nFinite-element twin (heat-sink sensor): measured coupling versus the 1D plane-source estimate. d = distance from the sensor to the")
print("nearest edge of the neighbouring footprint (half module width 43 mm + gap).")
for gap, tag in ((25, "fem_g25_h{h}_rb1_src1.csv"), (40, "fem_g40_h{h}_rb1_src1.csv"), (55, "fem_g55_h{h}_rb1_src1.csv")):
    d = 0.043 + gap / 1000
    for h in (1000, 8000):
        S = np.loadtxt(os.path.join(fem, tag.format(h=h)), delimiter=",")[:, :8]
        rise = S[ts, 4:].mean(1) / S[ts, :4].mean(1)
        sl = np.gradient(S, axis=0); slope = sl[ts, 4:].mean(1) / sl[ts, :4].mean(1)
        ok = np.all(rise <= c_rise(ts, d, a_al) + 1e-12) and np.all(slope <= c_slope(ts, d, a_al) + 1e-12)
        print(f" gap {gap} mm, h {h}: Fo(4 s) {a_al*4/d**2:.4f}")
        print("   rise  coupling FEM  " + " ".join(f"{x:.2e}" for x in rise) + "   1D " + " ".join(f"{x:.2e}" for x in c_rise(ts, d, a_al)))
        print("   slope coupling FEM  " + " ".join(f"{x:.2e}" for x in slope) + "   1D " + " ".join(f"{x:.2e}" for x in c_slope(ts, d, a_al)))
        print(f"   FEM below the 1D estimate at all listed times ({ts.tolist()} s): {ok}")

# real board: fit the single dimensionless parameter d^2/alpha on the early rise-coupling curve, check the shape
y = np.load(sys.argv[2])["Sat"]; T, P = y[:8], y[8:]
s = 7201                                             # Q_h switched on, Q_l off (E25)
tt = np.arange(1, 61)
cr = (T[1, s + tt] - T[1, s - 1]) / (T[0, s + tt] - T[0, s - 1])
grid = np.geomspace(0.1, 1000, 4000)                 # tau = d^2 / alpha in seconds
sse = [np.sum((c_rise(tt, np.sqrt(g), 1.0) - cr) ** 2) for g in grid]
g = grid[int(np.argmin(sse))]
pred = c_rise(tt, np.sqrt(g), 1.0)
print(f"\nReal board, Q_h step: best d^2/alpha = {g:.2f} s; rms misfit over 1-60 s {np.sqrt(np.min(sse)/len(tt)):.3f}")
for k in (1, 5, 10, 20, 30, 60):
    print(f"   t {k:2d} s: measured rise coupling {cr[k-1]:.3f}, 1D shape {pred[k-1]:.3f}")
print(f"   implied Fo of a 4 s slope window {4/g:.2f} (rig twin about 0.04); slope coupling at 4 s {c_slope(4, np.sqrt(g), 1.0):.3f}")
