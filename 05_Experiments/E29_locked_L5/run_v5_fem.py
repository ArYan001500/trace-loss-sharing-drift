"""Locked test L5F of method v5 (03_Protocols/METHOD_v5.md) on the finite-element twin.
As E28 fem_scenarios.py, plus: module gap drawn from {25, 40, 55} mm, 95 % interval coverage, method v5 only.
Drift scenarios on the finite-element twin (PDE Toolbox step responses, fem_twin.m). Independent of the
Cauer generator used to design methods v3/v4.
Per scenario: cooling h in {1000, 8000} W/m2K (static), baseline and current periods driven by the real mission
current of two different phases (loss ~ I^2, mean 450 W per module), module-1 loss scale U(0.85, 1.15) per
period, current-period interface resistance of module 1 x1.3 with probability 0.3. Four positions per module,
gain N(1, 0.02), lag U(0, 3) s, noise 0.3 K. Thermal coupling between the modules is included.
Sensor placement: 'sink' (heat-sink surface, below the grease) or 'base' (baseplate bottom, above the grease).
python fem_scenarios.py ROOT PLACEMENT SEED N"""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
FEM = os.path.join(here, "..", "E28_fem_twin")
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import lag
from load_mission import load
import method_v5

root, PLACE, SEED, N = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
off = 0 if PLACE == "sink" else 8
df = load(root, 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
for c in PH.values():
    df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
IMP = {}
for gp in (25, 40, 55):
    for h in (1000, 8000):
        for rb in ("1", "1.3"):
            for s in (1, 2):
                fn = f"fem_h{h}_rb{rb}_src{s}.csv" if gp == 40 else f"fem_g{gp}_h{h}_rb{rb}_src{s}.csv"
                S = np.loadtxt(os.path.join(FEM, fn), delimiter=",")[:, off:off + 8]
                IMP[(gp, h, rb, s)] = np.diff(S, axis=0)


def temps(I, gp, h, rb, ploss):
    p = (I / 1e3) ** 2
    P1, P2 = ploss * p, p
    sc = 450.0 / p.mean()
    out = np.zeros((len(I), 8))
    for s, P in ((1, P1), (2, P2)):
        x = sc * np.r_[P, P]                       # run the profile twice, keep the second pass
        for k in range(8):
            out[:, k] += np.convolve(x, IMP[(gp, h, rb, s)][:, k])[len(I):2 * len(I)]
    return 25.0 + out


rng = np.random.default_rng(SEED)
rows = []
for n in range(N):
    h = (1000, 8000)[rng.integers(2)]; gp = (25, 40, 55)[rng.integers(3)]
    pb, pn = rng.choice(list("UVW"), 2, replace=False)
    g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
    lb, ln = rng.uniform(0.85, 1.15), rng.uniform(0.85, 1.15)
    fault = rng.random() < 0.3
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy()
    Tb = temps(Ib, gp, h, "1", lb); Tn = temps(In, gp, h, "1.3" if fault else "1", ln)
    Tc_b, Tc_n = np.full(len(Ib), 25.0), np.full(len(In), 25.0)
    def pairs(T, Tc):
        return [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, k + 4 * j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(Tc))
                              for j in range(2)]) for k in range(4)]
    bp, np_ = pairs(Tb, Tc_b), pairs(Tn, Tc_n)
    truth = ln / lb
    r = method_v5.estimate_drift((Ib, Tc_b, bp), (In, Tc_n, np_))
    est = r.get("drift", np.nan); lo, hi = r.get("ci95", (np.nan, np.nan))
    rows.append((fault, h, gp, truth, est, r["abstain"], lo <= truth <= hi))
    print(f"{n:2d} gap {gp} h {h} fault {int(fault)} truth {truth:.3f} est {est:.3f} ({100*(est/truth-1):+.2f} %) "
          f"{'ABSTAIN ' + r.get('reason', '') if r['abstain'] else ''}", flush=True)
R = np.array(rows, dtype=float)
err = 100 * (R[:, 4] / R[:, 3] - 1); clean = R[:, 0] == 0; fl = ~clean; ab = R[:, 5] == 1; rep = clean & ~ab
D1 = np.mean(np.abs(err[rep]) <= 3); D2 = np.mean(R[rep, 6] == 1); D3 = np.mean(ab[fl] | (np.abs(err[fl]) <= 3)); D4 = ab[clean].mean()
print(f"\n{PLACE}: clean n {clean.sum()} reported {rep.sum()}; interface n {fl.sum()}")
for name, v, ok in (("D1 clean |err|<=3 %", D1, D1 >= 0.90), ("D2 coverage", D2, D2 >= 0.85), ("D3 interface", D3, D3 >= 0.80), ("D4 clean abstention", D4, D4 <= 0.10)):
    print(f"{name}: {v:.3f} -> {'PASS' if ok else 'FAIL'}")
print(f"clean |err| median {np.median(np.abs(err[rep])):.2f} %, max {np.max(np.abs(err[rep])):.2f} %")
