# E40 Post hoc audit of the coupling allowance on the locked L6 scenarios (not a locked test)

Method v6 and the E38 runners are used unchanged, with the locked seeds. The script checks that the locked estimates and abstention decisions are reproduced before it evaluates anything.

## Why
- The coupling allowance of v6 maps each interval end y through f_c^{-1}(y) with c_max = 0.25.
- This inversion is exact for a balanced baseline only. The observed drift is f_c(x_n) − f_c(x_b), and the inverse can fall on either side of the true log drift.
  - Example: rho_b = 1.15, rho_n = 1.10, c = 0.25, half-width 2e-5 → interval [0.95658, 0.97376] against a truth of 0.95652.
- A conservative alternative, under a bound |ln rho| <= L in both periods and a constant symmetric c <= c_max:
  - k = (1 − c^2)/(1 + c^2 + 2c cosh L);
  - envelope [min(A, A/k), max(B, B/k)] on the log interval [A, B].
  - With c = 0.25 and L = ln(1/0.85), the largest |ln rho| of the twin design, 1/k = 1.6737.

## Results (`interval_audit.py`, outputs `audit_{cauer,sink,base}.txt`; clean reported cases)

| run | n | coverage without allowance | with v6 allowance | with envelope | median half-width (none / v6 / envelope) |
|---|---|---|---|---|---|
| L6C Cauer | 81 | 0.975 | 1.000 | 1.000 | 1.57 / 4.95 / 4.98 % |
| L6F sink | 65 | 0.800 | 1.000 | 1.000 | 1.29 / 5.15 / 5.17 % |
| L6F base | 65 | 0.892 | 1.000 | 1.000 | 0.86 / 5.06 / 5.08 % |

Reading:
- The allowance is needed under coupling: without it, coverage falls to 0.80 below the grease.
- On these scenarios, the v6 allowance and the conservative envelope give the same coverage and nearly the same width.
- The envelope rests on an assumed bound L on the absolute loss ratio, which temperatures cannot reveal (Proposition 1). It is therefore a conditional bound, not a data-driven guarantee.

Command: `python interval_audit.py ROOT cauer | sink | base` (ROOT = extracted Zenodo record 20395811).
