"""Coupling seen by the guard statistics of TRACE in the finite-element twin with the corrected grease (E37);
same computation as E36 guard_window_coupling.py on the E28 twin.
For a loss step in module 1 at t = 0, the slope jump of Eq. (15) over a window [a, b) is the mean central
difference of the rise over samples a..b-1 (the pre-step slope is zero). The ratio neighbor/own of these
jumps is the coupling seen by that window. Windows: [1,4) and [5,9) for G1, [1,5) for G2 and the slope-jump drift.
python guard_window_coupling.py  (run from this folder)"""
import os
import numpy as np

FEM = os.path.dirname(os.path.abspath(__file__))


def jump(y, a, b):
    return np.gradient(y)[a:b].mean()


rows = []
for gap in (25, 40, 55):
    for h in (1000, 8000):
        fn = f"fem_g{gap}_h{h}_rb1_src1.csv"
        S = np.loadtxt(os.path.join(FEM, fn), delimiter=",")
        for place, own, nb in (("below grease", slice(0, 4), slice(4, 8)), ("above grease", slice(8, 12), slice(12, 16))):
            o, n = S[:, own].mean(1), S[:, nb].mean(1)
            rows.append((gap, h, place, *[jump(n, a, b) / jump(o, a, b) for a, b in ((1, 4), (1, 5), (5, 9))]))

print("gap  h      placement      [1,4)      [1,5)      [5,9)")
for r in rows:
    print("%3d %5d  %-13s %9.2e  %9.2e  %9.2e" % r)
print("maximum over all: [1,4) %.2e  [1,5) %.2e  [5,9) %.2e" % tuple(max(r[i] for r in rows) for i in (3, 4, 5)))
