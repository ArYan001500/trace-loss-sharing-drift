"""Cross-implementation check. Exports test cases to CSV, computes the transfer drift with the Python code
(src/tdrift.py), runs the independent MATLAB implementation (tdrift_m.m), and compares.
Cases: rig mission halves (U, V, W) and 20 development scenarios of the Cauer generator (seed 7, draws as L3).
python export_and_compare.py ROOT MATLAB_EXE"""
import sys, os, json, subprocess
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
import tdrift
from load_mission import load
from cauer_ds import simulate, lag

root, mexe = sys.argv[1], sys.argv[2]
cases = []
df = load(root, 0); n = len(df); h = n // 2
CIRC = {"U": "A", "V": "A", "W": "C"}
for ph in "UVW":
    Tc = (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
    ks = [k for k in "1234" if f"{ph}1T{k}" in df and f"{ph}2T{k}" in df]
    A1 = np.column_stack([df[f"{ph}1T{k}"] for k in ks]); A2 = np.column_stack([df[f"{ph}2T{k}"] for k in ks])
    cases.append((f"rig_{ph}", Tc[:h], A1[:h], A2[:h], Tc[h:], A1[h:], A2[h:]))
PAR = json.load(open(os.path.join(here, "..", "E16_datasheet_L3", "ds_params.json")))
d2 = df.copy(); PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
for c in PH.values(): d2[c] = d2[c].shift(-11)
d2 = d2.dropna(subset=list(PH.values())).reset_index(drop=True)
rng = np.random.default_rng(7)
for s in range(20):
    lev = ["0.5", "1.0", "2.0"][rng.integers(3)]; pb, pn = rng.choice(list("UVW"), 2, replace=False)
    cap = rng.uniform(0.95, 1.05); g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
    ab = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.8, 1.25), cap=cap, rb=1.0)
    an = dict(ploss=rng.uniform(0.85, 1.15), rh=rng.uniform(0.8, 1.25), cap=cap, rb=1.0 if rng.random() < 0.7 else rng.uniform(1.15, 1.4))
    out = []
    for p, a in ((pb, ab), (pn, an)):
        I = d2[PH[p]].to_numpy(); Tc = (0.5 * (d2[f"Tw_{CIRC[p]}-In"] + d2[f"Tw_{CIRC[p]}-Out"])).to_numpy()
        T, _, _ = simulate(I, Tc, PAR[lev], a)
        R = [np.column_stack([lag(Tc + g[k, j] * (T[:, j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(I)) for k in range(4)]) for j in range(2)]
        out.append((Tc, R[0], R[1]))
    cases.append((f"sim_{s:02d}", *out[0], *out[1]))
os.makedirs(os.path.join(here, "cases"), exist_ok=True)
py = {}
for name, Tcb, B1, B2, Tcn, N1, N2 in cases:
    for tag, arr in (("Tcb", Tcb[:, None]), ("B1", B1), ("B2", B2), ("Tcn", Tcn[:, None]), ("N1", N1), ("N2", N2)):
        np.savetxt(os.path.join(here, "cases", f"{name}_{tag}.csv"), arr, delimiter=",", fmt="%.10f")
    r = tdrift.estimate_drift(Tcb, [(f"T{k}", B1[:, k], B2[:, k]) for k in range(B1.shape[1])],
                              Tcn, [(f"T{k}", N1[:, k], N2[:, k]) for k in range(N1.shape[1])])
    # Python works on the CSV-rounded values for a like-for-like comparison
    rd = lambda t: np.loadtxt(os.path.join(here, "cases", f"{name}_{t}.csv"), delimiter=",", ndmin=2)
    Tcb, B1, B2, Tcn, N1, N2 = rd("Tcb")[:, 0], rd("B1"), rd("B2"), rd("Tcn")[:, 0], rd("N1"), rd("N2")
    r = tdrift.estimate_drift(Tcb, [(f"T{k}", B1[:, k], B2[:, k]) for k in range(B1.shape[1])],
                              Tcn, [(f"T{k}", N1[:, k], N2[:, k]) for k in range(N1.shape[1])])
    py[name] = (r["drift"], *r["ci95"])
names = [c[0] for c in cases]
script = "names = {" + ",".join(f"'{n}'" for n in names) + "}; fid = fopen('matlab_results.csv','w'); " \
         "for i = 1:numel(names); r = @(t) readmatrix(['cases/' names{i} '_' t '.csv']); " \
         "[d, ci] = tdrift_m(r('Tcb'), r('B1'), r('B2'), r('Tcn'), r('N1'), r('N2')); " \
         "fprintf(fid, '%s,%.12f,%.12f,%.12f\\n', names{i}, d, ci(1), ci(2)); end; fclose(fid);"
subprocess.run([mexe, "-batch", script], cwd=here, check=True)
mat = {l.split(",")[0]: tuple(map(float, l.strip().split(",")[1:])) for l in open(os.path.join(here, "matlab_results.csv"))}
worst = 0.0
for n_ in names:
    p, m = np.array(py[n_]), np.array(mat[n_])
    dev = np.max(np.abs(np.log(p) - np.log(m))); worst = max(worst, dev)
    print(f"{n_}: python {p[0]:.6f} [{p[1]:.6f}, {p[2]:.6f}]  matlab {m[0]:.6f} [{m[1]:.6f}, {m[2]:.6f}]  max |log diff| {dev:.2e}")
print(f"\nworst |log difference| over {len(names)} cases (point and both interval ends): {worst:.2e}")
