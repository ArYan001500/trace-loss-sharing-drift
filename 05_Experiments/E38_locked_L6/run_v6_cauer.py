"""Test of method v6 (03_Protocols/METHOD_v6.md) on the Cauer generator. Draw protocol of L3 and L5C
(E16 run_drift_ds.py, E29 run_v5_cauer.py), plus one draw per scenario after the others: with probability 0.5 the
module labels are swapped in both periods, so that an interface fault falls on the reference module (module 2);
the true drift is then inverted.
python run_v6_cauer.py ROOT SEED N OUT [LAGMAX]"""
import sys, os, json
import numpy as np, pandas as pd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import simulate, lag
from load_mission import load
import method_v6
from common import swap, summary

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
    side = 2 if rng.random() < 0.5 else 1
    Ib, Tcb = drive(pb)
    In, Tcn = drive(pn)
    Tb, _, Lb = simulate(Ib, Tcb, PAR[lev], ab)
    Tn, _, Ln = simulate(In, Tcn, PAR[lev], an)
    rb_, rn_ = readings(Ib, Tcb, Tb, g, tl, rng), readings(In, Tcn, Tn, g, tl, rng)
    td = truth(In, Ln) / truth(Ib, Lb)
    if side == 2:
        rb_, rn_, td = swap(rb_), swap(rn_), 1 / td
    r = method_v6.estimate_drift((Ib, Tcb, rb_), (In, Tcn, rn_))
    est = r.get("drift", np.nan)
    lo, hi = r.get("ci95", (np.nan, np.nan))
    rows.append(dict(s=s, lev=lev, pb=pb, pn=pn, cap=cap, rb=an["rb"], side=side, fault=an["rb"] > 1.0, truth=td,
                     est=est, err_pct=100 * (est / td - 1), covered=bool(lo <= td <= hi), abstain=r["abstain"],
                     code="|".join(r["code"]), reason=r.get("reason", ""), zx=r.get("zx", np.nan),
                     shape_z=r.get("shape_z", np.nan), v3=r.get("v3", np.nan), lo=lo, hi=hi))
    q = rows[-1]
    print(f"{s:2d} {lev} rb {an['rb']:.2f} side {side} truth {td:.3f} est {est:.3f} err {q['err_pct']:+.2f} % "
          f"zx {q['zx']:+.2f} shz {q['shape_z']:+.2f} {'ABSTAIN ' + r['reason'] if r['abstain'] else ''}", flush=True)
t = pd.DataFrame(rows)
t.to_csv(OUT, index=False)
summary(t, "Cauer")
