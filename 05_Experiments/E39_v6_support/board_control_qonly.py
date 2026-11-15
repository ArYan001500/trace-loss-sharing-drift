"""Negative control of the Fourier-number domain rule on real hardware dynamics.
The linear thermal model of the two-transistor buck board is identified from its measured powers and
temperatures (public dataset, saturation experiment, E25; held-out RMSE 0.46 / 0.65 K). On this compact board the
neighbour's heat reaches the sensor within seconds (Fo of the slope window not small), so the method is outside
its stated domain. Expected: compressed drift; either abstention or failure of coverage.
Variant: one position only (transistor thermocouples). Scenarios as L3 without interface faults,
loss of each transistor ~ I^2 of the real mission current of a random phase (mean 1.5 W), loss scale of Q_h
U(0.85, 1.15) per period, noise 0.1 K, no added sensor lag. Frozen method v6, unchanged.
python board_control.py SANZ_NPZ ROOT SEED N"""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
import method_v6
from load_mission import load

S = np.load(sys.argv[1])["Sat"]
T, P = S[:8], S[8:]
W = T[:, 1:] @ np.linalg.pinv(np.vstack([T[:, :-1], P[:, :-1]]))
A, B = W[:, :8], W[:, 8:]
df = load(sys.argv[2], 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
for c in PH.values():
    df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
SEED, N = int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(SEED)


def run(I, rho):
    p = (I / 1e3) ** 2; p = 1.5 * p / p.mean()
    Pin = np.zeros((6, len(p))); Pin[0] = rho * p; Pin[1] = p
    x = np.linalg.solve(np.eye(8) - A, B @ Pin.mean(1)); out = np.empty((len(p), 8))
    for k in range(len(p)):
        out[k] = x; x = A @ x + B @ Pin[:, k]
    return 25.0 + out + rng.normal(0, 0.1, out.shape)


rows = []
for n in range(N):
    pb, pn = rng.choice(list("UVW"), 2, replace=False)
    lb, ln = rng.uniform(0.85, 1.15), rng.uniform(0.85, 1.15)
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy()
    Yb, Yn = run(Ib, lb), run(In, ln)
    pairs = lambda Y: [("Q", Y[:, 0], Y[:, 1])]
    Tc_b, Tc_n = np.full(len(Ib), 25.0), np.full(len(In), 25.0)
    r = method_v6.estimate_drift((Ib, Tc_b, pairs(Yb)), (In, Tc_n, pairs(Yn)))
    truth = ln / lb
    est = r.get("drift", np.nan); lo, hi = r.get("ci95", (np.nan, np.nan))
    rows.append((truth, est, r["abstain"], lo <= truth <= hi, r.get("v3", np.nan)))
    print(f"{n:2d} truth {truth:.3f} est {est:.3f} ({100*(est/truth-1):+.1f} %) v3 {r.get('v3', np.nan):.3f} "
          f"ci [{lo:.3f}, {hi:.3f}] {'ABSTAIN ' + r.get('reason', '') if r['abstain'] else ''}", flush=True)
R = np.array(rows, dtype=float)
err = 100 * (R[:, 1] / R[:, 0] - 1); ab = R[:, 2] == 1; rep = ~ab
slope = np.polyfit(np.log(R[:, 0]), np.log(R[:, 1]), 1)[0]
print(f"\nn {N}: abstained {ab.mean():.3f}; reported: within 3 % {np.mean(np.abs(err[rep]) <= 3):.3f}, "
      f"coverage {np.mean(R[rep, 3] == 1):.3f}, median |err| {np.median(np.abs(err[rep])):.2f} %; "
      f"compression slope d(log est)/d(log truth) {slope:.3f} (1 = no compression)")
