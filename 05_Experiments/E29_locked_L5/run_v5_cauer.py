"""Locked test L5C of method v5 (03_Protocols/METHOD_v5.md), Cauer generator. Draw protocol identical to L3
(E16 run_drift_ds.py). python run_v5_cauer.py ROOT SEED N OUT [LAGMAX]"""
import sys, os, json
import numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import simulate, lag
from load_mission import load
import method_v5 as method_v4

root = sys.argv[1]
SEED, N, OUT = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
LAGMAX = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
PAR = json.load(open(os.path.join(here, "..", "E16_datasheet_L3", "ds_params.json")))
df = load(root, 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
CIRC = {"U": "A", "V": "A", "W": "C"}
for c in PH.values():
    df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)


def drive(ph):
    return df[PH[ph]].to_numpy(), (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()


def readings(I, Tc, T, g, tl, rng):
    return [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(I)) for j in range(2)])
            for k in range(4)]


def truth(I, L):
    m = I > np.median(I)
    return L[m, 0].mean() / L[m, 1].mean()


rng = np.random.default_rng(SEED)
rows = []
for s in range(N):
    lev = ["0.5", "1.0", "2.0"][rng.integers(3)]
    pb, pn = rng.choice(list("UVW"), 2, replace=False)
    cap = rng.uniform(0.95, 1.05)
    g = rng.normal(1, 0.02, (4, 2))
    tl = rng.uniform(0, LAGMAX, (4, 2))
    ab = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap, rb=1.0)
    an = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap,
              rb=1.0 if rng.random() < 0.7 else rng.uniform(1.15, 1.40))
    Ib, Tcb = drive(pb)
    In, Tcn = drive(pn)
    Tb, _, Lb = simulate(Ib, Tcb, PAR[lev], ab)
    Tn, _, Ln = simulate(In, Tcn, PAR[lev], an)
    rb_, rn_ = readings(Ib, Tcb, Tb, g, tl, rng), readings(In, Tcn, Tn, g, tl, rng)
    td = truth(In, Ln) / truth(Ib, Lb)
    r = method_v4.estimate_drift((Ib, Tcb, rb_), (In, Tcn, rn_))
    est = r.get("drift", np.nan)
    lo, hi = r.get("ci95", (np.nan, np.nan))
    rows.append(dict(s=s, lev=lev, pb=pb, pn=pn, cap=cap, rb=an["rb"], truth=td, est=est,
                     err_pct=100 * (est / td - 1), covered=lo <= td <= hi, abstain=r["abstain"],
                     reason=r.get("reason", ""), zx=r.get("zx", np.nan), shape_z=r.get("shape_z", np.nan),
                     v3=r.get("v3", np.nan)))
    print(f"{s:2d} {lev} rb {an['rb']:.2f} truth {td:.3f} est {est:.3f} err {rows[-1]['err_pct']:+.2f} % "
          f"zx {rows[-1]['zx']:+.2f} shz {rows[-1]['shape_z']:+.2f} {'ABSTAIN ' + r['reason'] if r['abstain'] else ''}", flush=True)
t = pd.DataFrame(rows)
t.to_csv(OUT, index=False)
clean = t[t.rb == 1.0]; rep = clean[~clean.abstain]; iface = t[t.rb > 1.0]
D1 = (rep.err_pct.abs() <= 3).mean(); D2 = rep.covered.mean()
D3 = (iface.abstain | (iface.err_pct.abs() <= 3)).mean(); D4 = clean.abstain.mean()
print(f"\nclean n={len(clean)} reported={len(rep)}; interface n={len(iface)}")
for name, val, ok in (("D1 clean |err|<=3 %", D1, D1 >= 0.90), ("D2 coverage, clean reported", D2, D2 >= 0.85),
                      ("D3 interface abstain or |err|<=3 %", D3, D3 >= 0.80), ("D4 clean abstention", D4, D4 <= 0.10)):
    print(f"{name}: {val:.3f} -> {'PASS' if ok else 'FAIL'}")
print(f"clean |err| median {rep.err_pct.abs().median():.2f} %, max {rep.err_pct.abs().max():.2f} %; "
      f"interface err median {iface.err_pct.median():+.1f} %")
