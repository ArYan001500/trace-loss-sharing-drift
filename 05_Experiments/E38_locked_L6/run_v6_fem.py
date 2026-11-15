"""Test of method v6 (03_Protocols/METHOD_v6.md) on the finite-element twin with the corrected grease resistance
(E37 step responses). Draw protocol of L5F (E29 run_v5_fem.py), plus one draw per scenario after the others: with
probability 0.5 the module labels are swapped in both periods, so that the interface fault falls on the reference
module (module 2); the true drift is then inverted.
Per scenario: cooling h in {1000, 8000} W/m2K, module gap in {25, 40, 55} mm, baseline and current periods driven
by the real mission current of two different phases (loss ~ I^2, mean 450 W per module), module-1 loss scale
U(0.85, 1.15) per period, current-period interface resistance of module 1 x1.3 with probability 0.3. Four positions
per module, gain N(1, 0.02), lag U(0, 3) s, noise 0.3 K. Thermal coupling between the modules is included.
Sensor placement: 'sink' (heat-sink surface, below the grease) or 'base' (baseplate bottom, above the grease).
python run_v6_fem.py ROOT PLACEMENT SEED N OUT"""
import sys, os
import numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
FEM = os.path.join(here, "..", "E37_fem_twin_grease")
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import lag
from load_mission import load
import method_v6
from common import swap, summary

root, PLACE, SEED, N, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
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
                S = np.loadtxt(os.path.join(FEM, f"fem_g{gp}_h{h}_rb{rb}_src{s}.csv"), delimiter=",")[:, off:off + 8]
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
    side = 2 if rng.random() < 0.5 else 1
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy()
    Tb = temps(Ib, gp, h, "1", lb); Tn = temps(In, gp, h, "1.3" if fault else "1", ln)
    Tc_b, Tc_n = np.full(len(Ib), 25.0), np.full(len(In), 25.0)

    def pairs(T, Tc):
        return [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, k + 4 * j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(Tc))
                              for j in range(2)]) for k in range(4)]
    bp, np_ = pairs(Tb, Tc_b), pairs(Tn, Tc_n)
    truth = ln / lb
    if side == 2:
        bp, np_, truth = swap(bp), swap(np_), 1 / truth
    r = method_v6.estimate_drift((Ib, Tc_b, bp), (In, Tc_n, np_))
    est = r.get("drift", np.nan); lo, hi = r.get("ci95", (np.nan, np.nan))
    rows.append(dict(n=n, gap=gp, h=h, pb=pb, pn=pn, side=side, fault=fault, truth=truth, est=est,
                     err_pct=100 * (est / truth - 1), covered=bool(lo <= truth <= hi), abstain=r["abstain"],
                     code="|".join(r["code"]), reason=r.get("reason", ""), zx=r.get("zx", np.nan),
                     shape_z=r.get("shape_z", np.nan), v3=r.get("v3", np.nan), lo=lo, hi=hi))
    q = rows[-1]
    print(f"{n:2d} gap {gp} h {h} fault {int(fault)} side {side} truth {truth:.3f} est {est:.3f} ({q['err_pct']:+.2f} %) "
          f"zx {q['zx']:+.2f} shz {q['shape_z']:+.2f} {'ABSTAIN ' + r.get('reason', '') if r['abstain'] else ''}", flush=True)
t = pd.DataFrame(rows)
t.to_csv(OUT, index=False)
summary(t, PLACE)
