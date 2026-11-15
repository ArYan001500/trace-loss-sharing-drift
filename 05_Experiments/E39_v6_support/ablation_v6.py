"""Ablation of the TRACE components (frozen method v6) on development scenarios.
Configurations, all with the same transfer-referenced point estimate and the v6 guard statistics:
  A  no guard, interval without coupling term
  B  shape guard only, two-sided (|z_s| > 2.2414)
  C  cross-band guard only, two-sided (|z_x| > 2.2414)
  D  both guards, one-sided as in method v5 (z_s < -1.96 or z_x > 1.96)
  E  both guards, two-sided (|z| > 2.2414)
  F  E + coupling interval                                   (= method v6, TRACE)
Generators and draws identical to the v6 runners (05_Experiments/E38_locked_L6): Cauer (development seed 7) and
finite-element twin with the corrected grease (E37; development seed 11; heat-sink and baseplate sensors); the
faulted module is drawn at random. Criteria as the locked tests: D1 clean within 3 %, D2 coverage, D3 interface
abstained or within 3 %, D4 clean abstention; last column the median interval half-width of clean reported cases.
python ablation_v6.py ROOT GEN SEED N   with GEN in {cauer, fem_sink, fem_base}"""
import sys, os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import simulate, lag
from load_mission import load
import method_v6

root, GEN, SEED, N = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
df = load(root, 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}; CIRC = {"U": "A", "V": "A", "W": "C"}
for c in PH.values(): df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
PAR = json.load(open(os.path.join(here, "..", "E16_datasheet_L3", "ds_params.json")))
FEM = os.path.join(here, "..", "E37_fem_twin_grease")
rng = np.random.default_rng(SEED)
coolant = lambda ph: (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
swap = lambda P: [(n, b, a) for n, a, b in P]
off = 8 if GEN == "fem_base" else 0
IMP = {}
if GEN.startswith("fem"):
    for gp in (25, 40, 55):
        for h in (1000, 8000):
            for rb in ("1", "1.3"):
                for s in (1, 2):
                    S = np.loadtxt(os.path.join(FEM, f"fem_g{gp}_h{h}_rb{rb}_src{s}.csv"), delimiter=",")
                    IMP[(gp, h, rb, s)] = np.diff(S[:, off:off + 8], axis=0)


def fem_temps(I, gp, h, rb, l):
    p = (I / 1e3) ** 2; sc = 450.0 / p.mean(); out = np.zeros((len(I), 8))
    for s, P in ((1, l * p), (2, p)):
        x = sc * np.r_[P, P]
        for k in range(8):
            out[:, k] += np.convolve(x, IMP[(gp, h, rb, s)][:, k])[len(I):2 * len(I)]
    return 25.0 + out


rows = []
for n in range(N):
    if GEN == "cauer":                                   # identical draw order to run_v6_cauer.py
        lev = ["0.5", "1.0", "2.0"][rng.integers(3)]; pb, pn = rng.choice(list("UVW"), 2, replace=False)
        cap = rng.uniform(0.95, 1.05); g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
        ab = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap, rb=1.0)
        an = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.80, 1.25), cap=cap,
                  rb=1.0 if rng.random() < 0.7 else rng.uniform(1.15, 1.40))
        side = 2 if rng.random() < 0.5 else 1
        Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy(); Tcb, Tcn = coolant(pb), coolant(pn)
        Tb, _, Lb = simulate(Ib, Tcb, PAR[lev], ab); Tn, _, Ln = simulate(In, Tcn, PAR[lev], an)
        rd = lambda I, Tc, T: [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(I)) for j in range(2)]) for k in range(4)]
        bp, np_ = rd(Ib, Tcb, Tb), rd(In, Tcn, Tn)
        tr = lambda I, L: L[I > np.median(I), 0].mean() / L[I > np.median(I), 1].mean()
        truth = tr(In, Ln) / tr(Ib, Lb); fault = an["rb"] > 1
    else:                                                # identical draw order to run_v6_fem.py
        h = (1000, 8000)[rng.integers(2)]; gp = (25, 40, 55)[rng.integers(3)]
        pb, pn = rng.choice(list("UVW"), 2, replace=False)
        g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
        lb, ln = rng.uniform(0.85, 1.15), rng.uniform(0.85, 1.15); fault = rng.random() < 0.3
        side = 2 if rng.random() < 0.5 else 1
        Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy()
        Tb = fem_temps(Ib, gp, h, "1", lb); Tn = fem_temps(In, gp, h, "1.3" if fault else "1", ln)
        Tcb, Tcn = np.full(len(Ib), 25.0), np.full(len(In), 25.0)
        rd = lambda T, Tc: [(f"T{k+1}", *[lag(Tc + g[k, j] * (T[:, k + 4 * j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(Tc)) for j in range(2)]) for k in range(4)]
        bp, np_ = rd(Tb, Tcb), rd(Tn, Tcn); truth = ln / lb
    if side == 2:
        bp, np_, truth = swap(bp), swap(np_), 1 / truth
    r = method_v6.estimate_drift((Ib, Tcb, bp), (In, Tcn, np_))
    if "drift" not in r:
        rows.append((fault, truth, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)); continue
    rows.append((fault, truth, r["drift"], *r["ci95_no_coupling"], *r["ci95"], r["zx"], r["shape_z"]))
R = np.array(rows, dtype=float)
fault, truth, d, lo0, hi0, lo1, hi1, zx, sz = R.T
ok = np.isfinite(d); err = 100 * (d / truth - 1); clean = (fault == 0) & ok; fl = (fault == 1) & ok
zg = method_v6.LOCK["z_guard"]
cfg = {"A": (lambda: np.zeros(len(d), bool), False),
       "B": (lambda: np.abs(sz) > zg, False),
       "C": (lambda: np.abs(zx) > zg, False),
       "D": (lambda: (sz < -1.96) | (zx > 1.96), False),
       "E": (lambda: (np.abs(sz) > zg) | (np.abs(zx) > zg), False),
       "F": (lambda: (np.abs(sz) > zg) | (np.abs(zx) > zg), True)}
print(f"{GEN}, seed {SEED}: clean {clean.sum()}, interface {fl.sum()}, no estimate {(~ok).sum()}")
print("cfg  D1 clean<=3%  D2 coverage  D3 interface  D4 clean abstain  median half-width %")
for name, (rule, coup) in cfg.items():
    ab = rule()
    lo, hi = (lo1, hi1) if coup else (lo0, hi0)
    cov = (lo <= truth) & (truth <= hi); rep = clean & ~ab
    print(f" {name}   {np.mean(np.abs(err[rep]) <= 3):10.3f} {np.mean(cov[rep]):11.3f} {np.mean(ab[fl] | (np.abs(err[fl]) <= 3)):13.3f} "
          f"{np.mean(ab[clean]):16.3f} {np.median(50 * (hi[rep] / lo[rep] - 1)):18.2f}")
