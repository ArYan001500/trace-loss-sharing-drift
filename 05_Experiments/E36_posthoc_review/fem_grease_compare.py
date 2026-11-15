"""Compare the finite-element twin with the grease area resistance of E28 (2.7e-4 m2K/W) and with the value
derived from the datasheet per half module over the full footprint (1.35e-4 m2K/W; run_fem_grease.m)."""
import numpy as np, os
here = os.path.dirname(os.path.abspath(__file__)); E28 = os.path.join(here, "..", "E28_fem_twin")
c = lambda S, o: S[-1, o + 4:o + 8].mean() / S[-1, o:o + 4].mean()
print("settled coupling at 900 s (below grease, above grease): E28 -> corrected")
for g in (25, 40, 55):
    for h in (1000, 8000):
        O = np.loadtxt(os.path.join(E28, f"fem_h{h}_rb1_src1.csv" if g == 40 else f"fem_g{g}_h{h}_rb1_src1.csv"), delimiter=",")
        N = np.loadtxt(os.path.join(here, f"fem_Rg_g{g}_h{h}_rb1_src1.csv"), delimiter=",")
        print(f"gap {g} h {h}: below {c(O, 0):.4f} -> {c(N, 0):.4f}   above {c(O, 8):.4f} -> {c(N, 8):.4f}")
print("interface resistance x1.3 under module 1 (gap 40): own response faulty/clean")
for h in (1000, 8000):
    for tag, F, C in (("E28", f"{E28}/fem_h{h}_rb1.3_src1.csv", f"{E28}/fem_h{h}_rb1_src1.csv"),
                      ("corrected", f"{here}/fem_Rg_g40_h{h}_rb1.3_src1.csv", f"{here}/fem_Rg_g40_h{h}_rb1_src1.csv")):
        F, C = np.loadtxt(F, delimiter=","), np.loadtxt(C, delimiter=",")
        r = lambda t, o: F[t, o:o + 4].mean() / C[t, o:o + 4].mean()
        print(f" h {h} {tag}: below 4 s {r(4, 0):.3f}, 900 s {r(900, 0):.3f} | above 4 s {r(4, 8):.3f}, 900 s {r(900, 8):.3f}")
