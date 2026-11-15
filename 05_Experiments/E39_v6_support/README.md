# E39 Supporting analyses with the frozen method v6 (development data and real recordings)

None of these runs is a locked test. Method v6 (`src/method_v6.py`, `03_Protocols/METHOD_v6.md`) is used unchanged.
ROOT = extracted Zenodo record 20395811; NPZ = `02_Data/working/sanz_1s/sanz_1s.npz` (board dataset, E25).

## Ablation (`ablation_v6.py`; Cauer seed 7, finite-element seed 11, 60 scenarios each)
Draws identical to the v6 runners (E38), including the random faulted module and the E37 twin.
Configuration F reproduces the E38 development runs exactly.
Each cell gives D1 / D2 / D3 / D4; the last value is the median interval half-width of clean reported cases.

| configuration | Cauer | FEM, sensor below grease | FEM, sensor above grease |
|---|---|---|---|
| A transfer drift only | 1.000 / 0.976 / **0.000** / 0.000; 1.55 % | 0.978 / **0.696** / 1.000 / 0.000; 1.29 % | 1.000 / **0.848** / **0.000** / 0.000; 0.85 % |
| B shape guard only, two-sided | 1.000 / 0.976 / 0.389 / 0.000; 1.55 % | 0.978 / 0.689 / 1.000 / 0.022; 1.31 % | 1.000 / 0.844 / 0.929 / 0.022; 0.85 % |
| C cross-band guard only, two-sided | 1.000 / 0.976 / 0.833 / 0.000; 1.55 % | 0.978 / 0.696 / 1.000 / 0.000; 1.29 % | 1.000 / 0.844 / 0.143 / 0.022; 0.85 % |
| D both guards, one-sided (v5 rule) | 1.000 / 0.976 / 0.278 / 0.000; 1.55 % | 0.978 / 0.689 / 1.000 / 0.022; 1.31 % | 1.000 / 0.860 / 0.357 / 0.065; 0.85 % |
| E both guards, two-sided | 1.000 / 0.976 / 0.833 / 0.000; 1.55 % | 0.978 / 0.689 / 1.000 / 0.022; 1.31 % | 1.000 / 0.841 / 0.929 / 0.043; 0.86 % |
| F E + coupling interval (= v6, TRACE) | 1.000 / 1.000 / 0.833 / 0.000; 5.10 % | 0.978 / 1.000 / 1.000 / 0.022; 5.20 % | 1.000 / 1.000 / 0.929 / 0.043; 4.36 % |

Reading:
- The two guards are complementary: the cross-band guard catches the Cauer faults (0.833 against 0.389 for the shape guard), and the shape guard catches the faults seen by the sensor above the grease (0.929 against 0.143).
- The one-sided rule of v5 (D) misses faults of module 2: D3 0.278 (Cauer) and 0.357 (above grease).
- With the sensor below the grease, the interface fault leaves the reported drift within 3 % even without guards (D3 1.000 for A).
- The coupling interval restores coverage (0.70 to 0.98 → 1.000). It costs width: the median half-width grows from 0.9 to 1.6 % to 4.4 to 5.2 %.

## Comparison with other estimators (`compare_baselines_v6.py`; seed 101, 40 scenarios per generator)
As E34, with v6, the E37 twin, and the faulted module drawn at random (labels swapped in losses, temperatures and calibration data).
Outputs: `cmp_cauer_output.txt`, `cmp_fem_output.txt`. Errors in %.

Cauer generator (27 clean, 13 interface: 10 on module 1, 3 on module 2):

| method | clean median | clean within 3 % | clean max | interface median, module 1 / module 2 |
|---|---|---|---|---|
| static ratio | 3.72 | 0.41 | 8.68 | +14.5 / −20.7 |
| wavelet ratio | 0.66 | 1.00 | 1.32 | +16.4 / −16.2 |
| two-temperature fit | 0.53 | 0.96 | 3.75 | +18.0 / −19.7 |
| slope jump (v3) | 1.36 | 0.89 | 3.47 | +5.9 / −10.1 |
| v6 | 0.41 | 1.00 | 1.99 | +15.1 / −16.8 (abstained 10/10 and 2/3) |
| calibrated inversion (labels) | 2.79 | 0.56 | 8.10 | +14.4 / −20.3 |

Finite-element twin, sensor below the grease (31 clean, 9 interface: 5 on module 1, 4 on module 2):

| method | clean median | clean within 3 % | clean max | interface median, module 1 / module 2 |
|---|---|---|---|---|
| static ratio | 1.02 | 0.74 | 8.80 | +1.2 / +1.9 |
| wavelet ratio | 0.37 | 1.00 | 1.55 | −0.6 / +0.9 |
| two-temperature fit | 0.35 | 1.00 | 2.05 | +0.2 / −0.1 |
| slope jump (v3) | 1.48 | 0.87 | 5.15 | −14.3 / +11.3 |
| v6 | 0.53 | 0.94 | 4.03 | +0.4 / +1.2 (abstained 5/5 and 4/4) |
| calibrated inversion (labels) | 0.56 | 0.97 | 4.20 | +0.3 / −1.0 |

v6 on reported cases only: clean within 3 % 1.00 (Cauer) and 0.97 (twin); clean abstention 0.00 and 0.06.

Reading:
- On the Cauer generator every other estimator reports the interface faults of either module as loss changes of 6 to 21 %, with no warning. v6 abstains on 12 of 13.
- On the twin with the sensor below the grease, the interface fault barely moves the wideband estimators (within 2 %). Only the slope jump is biased (−14 % and +11 %).

## Rig recording (`rig_v6_halves.py`, `rig_sensitivity.py`)
- First half against second half, per phase, both directions: all six cases reported, every interval contains 1, largest deviation 0.09 % (U, 1.0009).
- Intervals (first to second half): U [0.9849, 1.0186], V [0.9964, 1.0058], W [0.9906, 1.0073].
- Transfer drift over 27 analysis settings: 0.9955 to 1.0054 (U), 1.0002 to 1.0020 (V), 0.9933 to 1.0010 (W).
- Clock offset 10, 11 and 12 s: all reported, |z_x| <= 1.78, |z_s| <= 1.77.

## Board negative control (`board_control.py`, `board_control_qonly.py`; seed 41, 40 scenarios)
- Both variants: abstained 0.000; reported within 3 % 0.325, coverage 0.375, median |error| 3.85 %, compression slope 0.481.
- With two positions, v6 excludes the driver thermocouples because their stability statistic is non-finite (the slope-jump ratio is negative, checked on the first five scenarios), so both variants use the transistor thermocouples only.
- Method v5 abstained on all cases of the two-position variant because of a degenerate (NaN) guard value (E30). Under v6 this value no longer arises.
- Reading: outside the coupling condition, v6 reports about half of the true drift without warning, as v5 did on the one-position variant.

## Common loss waveform (`loss_mix_sensitivity.py`; copy of E36 on the E37 twin, gap 40 mm, h 1000 W/m²K, sensor below the grease, noise-free)
- Different mission phases in the two periods, current-proportional loss share up to 0.6: |error| <= 0.35 %.
- Current-period load scaled to 0.6 with a mixed loss law: errors −3.3 to −6.5 %; scaled to 0.3: −9.3 to −14.7 %.
- With a purely quadratic loss (m = 0), every case gives an error of 0.00 %.

## Worked example (`worked_example_v6.py`, `worked_example_v6_output.txt`)
Cauer development scenario 6 of seed 7 (L3 draw order), truth 1.0894.
- Position T1 is excluded by the stability test; T2, T3 and T4 are used.
- Random effects over the three: pooled log d = +0.08308, se = 0.00565, t(0.975, 2) = 4.303. Drift 1.0866 [1.0605, 1.1134].
- Guards: z_x −0.03, z_s −0.02.
- Coupling: f^-1(ln 1.1134) = +0.17930 gives the upper end 1.1964; the lower end stays at 1.0605.
- Error −0.25 %.
