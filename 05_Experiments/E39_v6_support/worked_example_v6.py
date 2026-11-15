"""Worked numeric example of TRACE (method v6) on one development scenario of the Cauer generator
(seed 7, scenario 6 of the L3 draw order). Prints every intermediate quantity so that the combination
step can be checked by hand.
python worked_example.py ROOT"""
import sys, os, json
import numpy as np
from scipy.stats import t as student_t
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import simulate, lag
from load_mission import load
import tdrift, method_v6

root = sys.argv[1]; TARGET = 6
df = load(root, 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}; CIRC = {"U": "A", "V": "A", "W": "C"}
for c in PH.values(): df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
PAR = json.load(open(os.path.join(here, "..", "E16_datasheet_L3", "ds_params.json")))
coolant = lambda ph: (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
rng = np.random.default_rng(7)
for n in range(TARGET + 1):
    lev = ["0.5", "1.0", "2.0"][rng.integers(3)]; pb, pn = rng.choice(list("UVW"), 2, replace=False)
    cap = rng.uniform(0.95, 1.05); g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
    ab = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap, rb=1.0)
    an = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap, rb=1.0 if rng.random() < 0.7 else rng.uniform(1.15, 1.40))
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy(); Tcb, Tcn = coolant(pb), coolant(pn)
    Tb, _, Lb = simulate(Ib, Tcb, PAR[lev], ab); Tn, _, Ln = simulate(In, Tcn, PAR[lev], an)
    rd = lambda I, Tc, T: [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(I)) for j in range(2)]) for k in range(4)]
    bp, np_ = rd(Ib, Tcb, Tb), rd(In, Tcn, Tn)
tr = lambda I, L: L[I > np.median(I), 0].mean() / L[I > np.median(I), 1].mean()
print(f"Scenario: case-capacity level {lev}, baseline phase {pb}, current phase {pn}")
print(f"  module-1 loss scale {ab['ploss']:.4f} -> {an['ploss']:.4f}; cooling factor {ab['rh']:.3f} -> {an['rh']:.3f}; capacity factor {cap:.3f}; interface {an['rb']:.2f}")
print(f"  sensor lags (s), module 1: {np.round(tl[:, 0], 2)}, module 2: {np.round(tl[:, 1], 2)}")
print(f"  truth: loss ratio baseline {tr(Ib, Lb):.4f}, now {tr(In, Ln):.4f}, drift {tr(In, Ln) / tr(Ib, Lb):.4f}\n")
r = method_v6.estimate_drift((Ib, Tcb, bp), (In, Tcn, np_))
print("Step 0, load steps and stability test (v6 uses one sensor set):", r["n_steps"], "steps; positions used", r["positions"], "; flags", r["flags"])
use = set(r["positions"])
s = tdrift.LOCK["skip"]
print("Step 1-3, per position: baseline map (6 Laguerre terms + direct term), scalar gain now, forward and reverse")
ys, vs = [], []
for (nm, b1, b2), (_, n1, n2) in [(a, b) for a, b in zip(bp, np_) if a[0] in use]:
    B1, B2, N1, N2 = b1 - Tcb, b2 - Tcb, n1 - Tcn, n2 - Tcn
    ib, i_n = np.arange(s, len(B1)), np.arange(s, len(N1))
    f = tdrift._one_dir(B2, B1, N2, N1, ib, i_n)["lev"]; r = tdrift._one_dir(B1, B2, N1, N2, ib, i_n)["lev"]
    p = tdrift.position(B1, B2, N1, N2)
    ys.append(p["log_d"]); vs.append(p["v"])
    print(f"  {nm}: forward gain {f:.5f}, reverse gain {r:.5f}, log d = (ln {f:.5f} - ln {r:.5f})/2 = {p['log_d']:+.5f} "
          f"(d = {np.exp(p['log_d']):.4f}), jackknife sd {np.sqrt(p['v']):.5f}")
y, v = np.array(ys), np.array(vs)
w = 1 / v; m0 = (w @ y) / w.sum(); Q = w @ (y - m0) ** 2; dfree = len(y) - 1
tau2 = max(0.0, (Q - dfree) / (w.sum() - (w ** 2).sum() / w.sum()))
ws = 1 / (v + tau2); mu = (ws @ y) / ws.sum(); q = (ws @ (y - mu) ** 2) / dfree
se = np.sqrt(max(q, 1.0) / ws.sum()); crit = student_t.ppf(0.975, dfree)
print(f"\nStep 4, random effects: Q = {Q:.3f} (df {dfree}), tau^2 = {tau2:.2e}, pooled log d = {mu:+.5f}, HKSJ factor q = {q:.3f} "
      f"(floored at 1), se = {se:.5f}, t(0.975, {dfree}) = {crit:.3f}")
print(f"  drift = exp({mu:+.5f}) = {np.exp(mu):.4f}; 95 % interval exp({mu:+.5f} -/+ {crit*se:.5f}) = [{np.exp(mu-crit*se):.4f}, {np.exp(mu+crit*se):.4f}]")
r = method_v6.estimate_drift((Ib, Tcb, bp), (In, Tcn, np_))   # recomputed: r was reused above
zg = method_v6.LOCK['z_guard']; c = method_v6.LOCK['c_max']
A, B = mu - crit * se, mu + crit * se
print(f"\nStep 5, guards: slope-jump drift {r['v3']:.4f}; z_x = {r['zx']:+.2f}, shape z = {r['shape_z']:+.2f} (abstain if either |z| > {zg:.3f})")
print(f"Step 6, coupling: ln A = {A:+.5f}, ln B = {B:+.5f}; f^-1(y) = ln((e^y - {c})/(1 - {c} e^y)); f^-1(A) = {method_v6._uncouple(A, c):+.5f}, f^-1(B) = {method_v6._uncouple(B, c):+.5f}; widened outward only")
print(f"\nTRACE output: drift {r['drift']:.4f}, 95 % interval [{r['ci95'][0]:.4f}, {r['ci95'][1]:.4f}], abstain {r['abstain']}; "
      f"error against truth {100 * (r['drift'] / (tr(In, Ln) / tr(Ib, Lb)) - 1):+.2f} %")
