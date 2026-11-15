# Method v6: two-sided guards, one sensor set, and an end-point coupling interval

Status: FROZEN on 2026-11-14, before locked tests L6C and L6F were run.
- Code: `src/method_v6.py`, `src/tdrift.py`, `src/method_v3.py`, `src/fingerprint_est.py`.
- Generators: `src/cauer_ds.py` with `05_Experiments/E16_datasheet_L3/ds_params.json`; finite-element step responses `05_Experiments/E37_fem_twin_grease/fem_g*_h*_rb*_src*.csv` from `fem_twin_grease.m` and `run_all.m`.
- Runners: `05_Experiments/E38_locked_L6/run_v6_cauer.py`, `run_v6_fem.py`, `common.py`.
- Hashes: `03_Protocols/FREEZE_v6_SHA256.txt`.

## Why v6
- A post hoc check of the locked L5 scenarios (`05_Experiments/E36_posthoc_review/label_swap.py`) swapped the module labels.
  - An interface fault then falls on the reference module (module 2).
  - Method v5 abstained on none of these faults: 0/20 (Cauer), 0/18 (sink), 0/18 (base).
  - Cause: both v5 guards were one-sided, and a fault of module 2 moves both statistics the other way.
- The v5 implementation differed from its description in three points:
  - an unstable position was excluded only from the slope-jump branch;
  - the guards used the raw temperature, while the point estimate used the coolant-referenced rise;
  - a non-finite z_x did not lead to abstention.
- The grease area resistance of the finite-element twin of E28 was twice the value implied by the datasheet.
  - E37 rebuilds all 24 step responses with R'' = 0.0225 x 0.086 x 0.070 = 1.35e-4 m²K/W.
  - Settled coupling, sensor below the grease, gap 25 mm, h 1000 W/m²K: 0.245 (E28: 0.242).
  - Settled coupling, sensor above the grease, same case: 0.200 (E28: 0.167).
  - c_max = 0.25 still bounds both.

## Estimator (changes from v5)
- V1, two-sided guards:
  - abstain if |z_s| > 2.2414 or |z_x| > 2.2414;
  - this is one family of two two-sided tests at 5 %, each at 2.5 % (Bonferroni).
- V2, one sensor set: a position with |z_stab| > 2 in either period is excluded from the guards and from the point estimate.
- V3, coolant-referenced rise: the guards use T − Tc, as the point estimate does.
- V4, invalid statistics: any non-finite statistic leads to abstention with a reason code:
  - `insufficient_excitation`, `no_stable_position`, `invalid_statistic`;
  - guard codes: `shape_change`, `band_disagreement`.
- V5, coupling interval by end-point mapping:
  - settled coupling c compresses the wideband log drift as f_c(l) = ln((e^l + c)/(1 + c e^l));
  - each end point y of the uncorrected interval is mapped with f_{c_max}^{-1}(y) = ln((e^y − c)/(1 − c e^y));
  - the interval is widened only outward: [min(A, f^{-1}(A)), max(B, f^{-1}(B))], with c_max = 0.25;
  - an end point outside the domain |y| < ln(1/c_max) is unbounded.
- Unchanged: point estimate (tdrift.py), slope-jump windows, load-step detection, stability test, bootstrap, and random-effects pooling.

## Development record (not locked evidence)
- Development seeds only: Cauer 7, finite-element twin 11; 60 scenarios per run.
  - Cauer: D1 1.000, D2 1.000, D3 0.833, D4 0.000.
    - Module-1 faults: 5 abstained, 1 reported beyond 3 %.
    - Module-2 faults: 10 abstained, 2 reported beyond 3 %.
  - Sink: D1 0.978, D2 1.000, D3 1.000, D4 0.022.
  - Base: D1 1.000, D2 1.000, D3 0.929, D4 0.043.
- Missed faults had guard statistics between 1.74 and 2.12 in magnitude.
- No threshold or setting was changed after these runs.
- Outputs: `05_Experiments/E38_locked_L6/dev_*`.

## Locked tests (criteria fixed now; each runs once)
- L6C, Cauer generator: `python run_v6_cauer.py ROOT 2041 100 locked_L6C_table.csv 3.0`.
- L6F, finite-element twin with the corrected grease, both placements:
  - `python run_v6_fem.py ROOT sink 2042 100 locked_L6F_sink_table.csv`
  - `python run_v6_fem.py ROOT base 2042 100 locked_L6F_base_table.csv`
- Draws as L5C and L5F, plus one draw per scenario: with probability 0.5 the module labels are swapped in both periods, so the interface fault falls on module 2 and the true drift is inverted.
- Criteria, for each of the three runs (unchanged from v5):
  - D1: clean reported cases with |error| <= 3 %, at least 0.90;
  - D2: 95 % interval coverage, clean reported, at least 0.85;
  - D3: interface cases abstained or |error| <= 3 %, at least 0.80;
  - D4: clean abstention, at most 0.10.
- v6 passes only if all three runs pass all four criteria.
- Reported in addition, without a criterion: the D3 counts per faulted module (abstained, reported within 3 %, reported beyond 3 %).
- Results are appended to `05_Experiments/E38_locked_L6/README.md` without changing code or protocol.

## Not claimed
- Real-hardware accuracy with a known loss-ratio reference.
- Separation of a gain-only path change from a loss change: both scale the transfer map alike and cannot be told apart.
- Guard calibration beyond the nominal level: the 5 % and 95 % levels are nominal.
- Coupling above 0.25, or a coupling bound independent of the twin used for the locked test.
