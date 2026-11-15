# Method v5: transfer-referenced loss-sharing drift with a family-wise guard and a coupling-aware interval

Status: FROZEN on 2026-11-13, before locked tests L5C and L5F were run.
- Code: `src/method_v5.py`, `src/tdrift.py`, `src/method_v3.py`, `src/fingerprint_est.py`.
- Generators: `src/cauer_ds.py` with `05_Experiments/E16_datasheet_L3/ds_params.json`; finite-element step responses `05_Experiments/E28_fem_twin/fem_*.csv` from `fem_twin.m`.
- Runners: `05_Experiments/E29_locked_L5/run_v5_cauer.py`, `run_v5_fem.py`.
- Hashes: `03_Protocols/FREEZE_v5_SHA256.txt`.

## Why v5
- Method v4 failed locked test L4 on D4 (clean abstention 0.116 against a maximum of 0.10).
  - Its two guards were each tested at 5 %, so the family-wise false-alarm rate was about twice the intended 5 %.
  - v5 tests each guard at 2.5 % (Bonferroni, |z| threshold 1.96).
  - This follows from the test design. It was not tuned on L4, whose data were not reused.
- A real two-device board (E25) and two independent twins of the rig (E27 finite volume, E28 finite element) show that thermal coupling between modules compresses the wideband drift towards 1.
  - The coupling is zero within 10 s on the rig layout.
  - Settled coupling reaches 0.24 for a module gap of 25 mm or more and cooling h of 1000 W/m²K or more.
  - v5 extends the interval on the side away from 1 by |log d| x 2 c_max / (1 − c_max), with c_max = 0.25.
  - The point estimate is unchanged.

## Estimator
- As method v4 (METHOD_v4.md sections 3 and 4), with three changes:
  - G1: abstain if the v3 shape index z < −1.96;
  - G2: abstain if z_x > 1.96;
  - the coupling extension of the interval described above.

## Development record (not locked evidence)
- Cauer generator, seeds 7, 8 and 9 (post hoc with the v5 thresholds): D1 1.000 / 1.000 / 1.000; D3 0.917 / 0.857 / 1.000; D4 0.021 / 0.000 / 0.026.
- Finite-element twin, seed 11, 60 per placement:
  - before the coupling term, sink coverage was 0.767;
  - after it: sink D1 1.000, D2 1.000, D3 1.000, D4 0.023; base D1 1.000, D2 1.000, D3 1.000, D4 0.000.

## Locked tests (criteria fixed now; each runs once)
- L5C, Cauer generator: `python run_v5_cauer.py ROOT 2032 60 locked_L5C_table.csv 3.0`. Draws identical to L3/L4.
- L5F, finite-element twin, both placements:
  - `python run_v5_fem.py ROOT sink 2031 60`
  - `python run_v5_fem.py ROOT base 2031 60`
  - Draws: module gap {25, 40, 55} mm, h {1000, 8000} W/m²K, real mission currents of two phases, module-1 loss scale U(0.85, 1.15) per period, interface x1.3 with probability 0.3, gain N(1, 0.02), lag U(0, 3) s, noise 0.3 K.
- Criteria, for each of the three runs:
  - D1: clean reported cases with |error| <= 3 %, at least 0.90;
  - D2: 95 % interval coverage, clean reported, at least 0.85;
  - D3: interface cases abstained or |error| <= 3 %, at least 0.80;
  - D4: clean abstention, at most 0.10.
- v5 passes only if all three runs pass all four criteria.
- Results are appended to `05_Experiments/E29_locked_L5/README.md` without changing code or protocol.

## Not claimed
- Real-hardware accuracy with a known loss-ratio reference.
- Geometry of the rig beyond figure-based estimates.
- Coupling above 0.25.
