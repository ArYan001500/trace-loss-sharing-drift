"""Real hardware (public buck-converter dataset, Sanz-Alcaine et al., TPEL 2024, saturation experiment):
two MOSFETs (high side Q_h, low side Q_l) on one board, each with a thermocouple, measured power per device.
Question: is the thermal cross-coupling between the two devices immediate (frequency-flat) or delayed?
At every power step (block boundary) the step sizes dPh, dPl are known. Fit, over all steps,
  fast:  jump of dT/dt at each sensor = a_h dPh + a_l dPl  (slope 5..60 s after minus 120..5 s before)
  slow:  change of the settled rise (last 10 min of the block) = A_h dPh + A_l dPl
Coupling coefficient of sensor Q_h = (term of the other device) / (term of its own device), fast and slow.
Then the loss-step ratio dPh/dPl is estimated from temperatures only, assuming identical devices and the
reciprocal, symmetric coupling c: the raw ratio uses jumps or settled changes directly; the decoupled ratio
solves r_h = a (dPh + c dPl), r_l = a (c dPh + dPl). Input: 02_Data/working/sanz_1s/sanz_1s.npz (1 s means)."""
import sys
import numpy as np

d = np.load(sys.argv[1])
y = d["Sat"]
T, P = y[:8], y[8:]
Th, Tl, Ph, Pl = T[0], T[1], P[0], P[1]
tot = Ph + Pl
dp = np.abs(np.diff(Ph)) + np.abs(np.diff(Pl))
cand = np.where(dp > 0.15)[0]
steps = []
for c in cand:
    if not steps or c - steps[-1] > 1200:
        steps.append(c)
rows = []
for s in steps:
    if s < 1300 or s + 3600 > len(Ph):
        continue
    nxt = [t for t in steps if t > s]
    end = (nxt[0] if nxt else len(Ph)) - 1
    if end - s < 3000:
        continue
    before = slice(s - 600, s)
    after = slice(end - 600, end)
    dPh, dPl = Ph[s + 600:s + 1200].mean() - Ph[before].mean(), Pl[s + 600:s + 1200].mean() - Pl[before].mean()
    sl = lambda x, a, b: np.polyfit(np.arange(a, b), x[s + a:s + b], 1)[0]
    jh, jl = sl(Th, 5, 60) - sl(Th, -120, -5), sl(Tl, 5, 60) - sl(Tl, -120, -5)
    Sh, Sl = Th[after].mean() - Th[before].mean(), Tl[after].mean() - Tl[before].mean()
    rows.append((s, dPh, dPl, jh, jl, Sh, Sl))
R = np.array(rows)
print(f"steps used: {len(R)} (block length >= 3000 s)")
X = R[:, 1:3]
for name, cols in (("fast (slope jump, K/s per W)", (3, 4)), ("slow (settled rise, K per W)", (5, 6))):
    ah = np.linalg.lstsq(X, R[:, cols[0]], rcond=None)[0]
    al = np.linalg.lstsq(X, R[:, cols[1]], rcond=None)[0]
    print(f"{name}: sensor Q_h own {ah[0]:.4g} other {ah[1]:.4g} -> coupling {ah[1]/ah[0]:.3f}; "
          f"sensor Q_l own {al[1]:.4g} other {al[0]:.4g} -> coupling {al[0]/al[1]:.3f}")
# ratio recovery on steps where both devices step in the same direction
m = (R[:, 1] * R[:, 2] > 0) & (np.abs(R[:, 2]) > 0.05)
true = R[m, 1] / R[m, 2]
if not m.any():
    print("ratio recovery: no step with both powers rising or falling; not evaluated"); sys.exit(0)
for name, (i, j) in (("fast", (3, 4)), ("slow", (5, 6))):
    raw = R[m, i] / R[m, j]
    ah = np.linalg.lstsq(X, R[:, i], rcond=None)[0]; al = np.linalg.lstsq(X, R[:, j], rcond=None)[0]
    c = 0.5 * (ah[1] / ah[0] + al[0] / al[1])
    dec = (raw - c) / (1 - c * raw)
    e_raw = 100 * (raw / true - 1); e_dec = 100 * (dec / true - 1)
    print(f"{name}: steps with both powers rising or falling: {m.sum()}; true ratio {true.min():.2f}..{true.max():.2f}; "
          f"raw median |err| {np.median(np.abs(e_raw)):.1f} %; decoupled (c = {c:.3f}, needs calibration) median |err| {np.median(np.abs(e_dec)):.1f} %")
    for t, r, dd in zip(true, raw, dec):
        print(f"   true {t:6.3f}  raw {r:6.3f}  decoupled {dd:6.3f}")
