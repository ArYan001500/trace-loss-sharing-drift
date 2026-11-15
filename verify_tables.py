"""Compares a re-computed case table with the stored reference table of a locked test, row by row:
scenario conditions, true drift, estimate, interval ends, coverage, abstention decision, reason codes, and the
two guard statistics; then recomputes criteria D1 to D4 from both tables.
python verify_tables.py NEW.csv REFERENCE.csv"""
import sys
import numpy as np, pandas as pd

a, b = pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2])
ok = len(a) == len(b) and list(a.columns) == list(b.columns)
num = ["truth", "est", "lo", "hi", "zx", "shape_z"]
same = ["fault", "side", "abstain", "covered", "code"] + [c for c in ("gap", "h", "pb", "pn", "lev") if c in a]
dev = {c: float(np.nanmax(np.abs(a[c] - b[c]))) for c in num}
eq = {c: bool((a[c].fillna("") == b[c].fillna("")).all()) for c in same}
ok = ok and all(v < 1e-9 for v in dev.values()) and all(eq.values())


def crit(t):
    clean, fault = ~t.fault, t.fault
    rep = clean & ~t.abstain
    e = t.err_pct.abs()
    return (round(float((e[rep] <= 3).mean()), 3), round(float(t.covered[rep].mean()), 3),
            round(float((t.abstain[fault] | (e[fault] <= 3)).mean()), 3), round(float(t.abstain[clean].mean()), 3))


ok = ok and crit(a) == crit(b)
print(f"{sys.argv[1]}: {len(a)} rows; identical fields {all(eq.values())}; "
      f"max |difference| {max(dev.values()):.1e}; D1-D4 {crit(a)} (reference {crit(b)}) -> {'MATCH' if ok else 'MISMATCH'}")
sys.exit(0 if ok else 1)
