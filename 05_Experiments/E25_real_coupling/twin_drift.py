"""Drift estimators under the REAL thermal coupling of a two-device board.
1. Identify a discrete linear thermal model T[k+1] = A T[k] + B P[k] (8 thermocouples, 6 measured power
   sources) from the saturation experiment of the public buck dataset (1 s means), as in the source paper.
   Check it open loop on the last 22 h of the same experiment, refitted without them.
2. Drive Q_h with rho * p(t) and Q_l with p(t), where p(t) follows the squared mission phase current of the
   traction rig, scaled to 1.5 W mean. Baseline rho_b = 1, current rho_n in {0.8, 0.9, 1.1, 1.25}.
   Sensors: the Q_h and Q_l thermocouples, white noise 0.1 K.
3. Estimators of rho_n / rho_b: transfer-referenced drift (wideband) and the slope-jump ratio drift
   (jump of dT/dt over samples 1-4 after each load step), both baseline-referenced."""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
import tdrift
from fingerprint_est import load_steps, jumps
from load_mission import load

d = np.load(sys.argv[1])
S, L = d["Sat"], d["Lin"]
T, P = S[:8], S[8:]
Z = np.vstack([T[:, :-1], P[:, :-1]])
W = T[:, 1:] @ np.linalg.pinv(Z)
A, B = W[:, :8], W[:, 8:]
print("max |eig A| =", np.abs(np.linalg.eigvals(A)).max().round(6))


def sim(Pin, x0):
    x = x0.copy(); out = np.empty((Pin.shape[1], 8))
    for k in range(Pin.shape[1]):
        out[k] = x; x = A @ x + B @ Pin[:, k]
    return out


# held-out check: refit on the first 50 h of the saturation experiment, simulate the last 22 h open loop
cut = 50 * 3600
W2 = T[:, 1:cut] @ np.linalg.pinv(np.vstack([T[:, :cut - 1], P[:, :cut - 1]]))
A2, B2 = W2[:, :8], W2[:, 8:]
x = T[:, cut].copy(); e = []
for k in range(cut, T.shape[1] - 1):
    e.append(x[:2] - T[:2, k]); x = A2 @ x + B2 @ P[:, k]
e = np.array(e)
print(f"held-out open-loop check (fit 0-50 h, simulate 50-73 h): RMSE Q_h {np.sqrt(np.mean(e[:, 0]**2)):.2f} K, "
      f"Q_l {np.sqrt(np.mean(e[:, 1]**2)):.2f} K; rise up to {T[0, cut:].max():.0f} / {T[1, cut:].max():.0f} K")
# single-source step responses of the identified model: coupling versus time
st = sim(np.vstack([np.ones(4000), np.zeros((5, 4000))]), np.zeros(8))
print("identified model, Q_h step: neighbour/own rise at 5, 30, 300, 3000 s:",
      [round(st[k, 1] / st[k, 0], 3) for k in (5, 30, 300, 3000)])

df = load(sys.argv[2], 0)
I = df["Ia_ef"].shift(-11).bfill().ffill().to_numpy()
p = (I / 1e3) ** 2; p = 1.5 * p / p.mean()
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
steps = load_steps(I, 0.15)


def run(rho):
    Pin = np.zeros((6, len(p))); Pin[0] = rho * p; Pin[1] = p
    # settle with the mean load first, then record
    x0 = np.linalg.solve(np.eye(8) - A, B @ Pin.mean(1))
    y = sim(Pin, x0)
    return y[:, 0] + rng.normal(0, 0.1, len(p)), y[:, 1] + rng.normal(0, 0.1, len(p))


def jr(a, b):
    j1, j2 = jumps(a, steps), jumps(b, steps); return (j1 @ j2) / (j2 @ j2)


zero = np.zeros(len(p))
for rep in range(3):
    b1, b2 = run(1.0)
    for rn in (0.8, 0.9, 1.1, 1.25):
        n1, n2 = run(rn)
        pos = tdrift.position(b1, b2, n1, n2)
        td = np.exp(pos["log_d"]); fj = jr(n1, n2) / jr(b1, b2)
        print(f"rep {rep} truth {rn:.3f}: transfer drift {td:.3f} ({100*(td/rn-1):+.1f} %), "
              f"slope-jump drift {fj:.3f} ({100*(fj/rn-1):+.1f} %)")
