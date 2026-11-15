"""Post hoc audit of the coupling allowance of method v6 on the locked scenarios of L6 (not a new locked test).
The frozen runners of E38 are executed unchanged with their locked seeds; every call of the frozen estimator is
recorded. For the clean reported cases, three intervals are compared:
  base      the transfer-referenced interval without coupling allowance (ci95_no_coupling)
  v6        the interval of method v6 (end-point mapping with c_max = 0.25)
  envelope  the derivative envelope [min(A, A/k), max(B, B/k)] in logs, k = (1 - c^2)/(1 + c^2 + 2c cosh L),
            c = 0.25, L = ln(1/0.85), the largest |ln rho| of the twin design (assumption, not enforced by data)
python interval_audit.py ROOT RUN   with RUN in {cauer, sink, base}"""
import sys, os, runpy
import numpy as np, pandas as pd

here = os.path.dirname(os.path.abspath(__file__))
E38 = os.path.join(here, "..", "E38_locked_L6")
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
import method_v6

root, run = sys.argv[1], sys.argv[2]
orig = method_v6.estimate_drift
log = []


def rec(base, now):
    r = orig(base, now)
    log.append(r)
    return r


method_v6.estimate_drift = rec
tmp = os.path.join(here, f"_tmp_{run}.csv")
if run == "cauer":
    sys.argv = ["run_v6_cauer.py", root, "2041", "100", tmp, "3.0"]
    script = "run_v6_cauer.py"
else:
    sys.argv = ["run_v6_fem.py", root, run, "2042", "100", tmp]
    script = "run_v6_fem.py"
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    runpy.run_path(os.path.join(E38, script), run_name="__main__")
t = pd.read_csv(tmp); os.remove(tmp)
ref = pd.read_csv(os.path.join(E38, f"locked_L6{'C' if run == 'cauer' else 'F_' + run}_table.csv"))
assert np.allclose(t.est, ref.est, equal_nan=True) and (t.abstain == ref.abstain).all(), "locked run not reproduced"

c, L = 0.25, np.log(1 / 0.85)
k = (1 - c * c) / (1 + c * c + 2 * c * np.cosh(L))
rows = []
for i, r in enumerate(log):
    if t.fault[i] or r["abstain"]:
        continue
    A, B = np.log(r["ci95_no_coupling"])
    env = (np.exp(min(A, A / k)), np.exp(max(B, B / k)))
    rows.append(dict(truth=t.truth[i], base=r["ci95_no_coupling"], v6=r["ci95"], env=env))
print(f"{run}: clean reported n = {len(rows)}; envelope factor 1/k = {1 / k:.4f} (c = {c}, L = ln(1/0.85))")
for name in ("base", "v6", "env"):
    cov = np.mean([q[name][0] <= q["truth"] <= q[name][1] for q in rows])
    hw = np.median([50 * (q[name][1] / q[name][0] - 1) for q in rows])
    print(f"  {name:4s}: coverage {cov:.3f}   median half-width {hw:.2f} %")
