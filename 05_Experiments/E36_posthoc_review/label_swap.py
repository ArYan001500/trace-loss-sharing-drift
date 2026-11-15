"""Post hoc check (not a locked test): the locked scenarios of L5C and L5F re-run with the module labels swapped.
The frozen runners of E29 are executed unchanged with their locked seeds. Every call of the frozen estimator
is evaluated twice on identical data: once as in the locked run, and once with modules 1 and 2 exchanged in both
periods. An interface fault of module 1 then appears as a fault of the reference module, and the true drift is
the inverse. Method v5 is not modified.
python label_swap.py DATA_ROOT RUN   with RUN in {cauer, sink, base}"""
import sys, os, runpy
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
E29 = os.path.join(here, "..", "E29_locked_L5")
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
import method_v5

root, run = sys.argv[1], sys.argv[2]
orig = method_v5.estimate_drift
log = []


def both(base, now):
    r = orig(base, now)
    sw = lambda per: (per[0], per[1], [(nm, b, a) for nm, a, b in per[2]])
    s = orig(sw(base), sw(now))
    log.append((r, s))
    return r


method_v5.estimate_drift = both
if run == "cauer":
    sys.argv = ["run_v5_cauer.py", root, "2032", "60", os.path.join(here, "_tmp_cauer.csv")]
    runpy.run_path(os.path.join(E29, "run_v5_cauer.py"), run_name="__main__")
    import pandas as pd
    t = pd.read_csv(os.path.join(here, "_tmp_cauer.csv")); os.remove(os.path.join(here, "_tmp_cauer.csv"))
    truth, fault = t.truth.to_numpy(), (t.rb > 1).to_numpy()
else:
    import io, contextlib, re
    buf = io.StringIO()
    sys.argv = ["run_v5_fem.py", root, run, "2031", "60"]
    with contextlib.redirect_stdout(buf):
        runpy.run_path(os.path.join(E29, "run_v5_fem.py"), run_name="__main__")
    lines = [l for l in buf.getvalue().splitlines() if re.match(r"\s*\d+ gap", l)]
    fault = np.array([int(re.search(r"fault (\d)", l).group(1)) == 1 for l in lines])
    truth = np.array([float(re.search(r"truth ([\d.]+)", l).group(1)) for l in lines])
    print(buf.getvalue().split("\n\n")[-1])


def summary(tag, res, tr):
    est = np.array([r.get("drift", np.nan) for r in res]); ab = np.array([r["abstain"] for r in res])
    lo = np.array([r.get("ci95", (np.nan, np.nan))[0] for r in res]); hi = np.array([r.get("ci95", (np.nan, np.nan))[1] for r in res])
    err = 100 * (est / tr - 1); clean = ~fault; rep = clean & ~ab
    D1 = np.mean(np.abs(err[rep]) <= 3); D2 = np.mean((lo[rep] <= tr[rep]) & (tr[rep] <= hi[rep]))
    D3 = np.mean(ab[fault] | (np.abs(err[fault]) <= 3)); D4 = ab[clean].mean()
    big = int(np.sum(~ab[fault] & (np.abs(err[fault]) > 3)))
    print(f"{tag}: D1 {D1:.3f}  D2 {D2:.3f}  D3 {D3:.3f}  D4 {D4:.3f}  | faults abstained {int(ab[fault].sum())}/{fault.sum()}, "
          f"reported beyond 3 %: {big}; clean abstained {int(ab[clean].sum())}/{clean.sum()}")
    if big:
        print("   errors of faults reported beyond 3 %:", np.round(err[fault & ~ab & (np.abs(err) > 3)], 1).tolist())


summary(f"{run} as locked ", [r for r, s in log], truth)
summary(f"{run} labels swapped", [s for r, s in log], 1 / truth)
