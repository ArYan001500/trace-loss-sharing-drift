# E25 Real thermal coupling between two devices (public buck-converter dataset)

## Data
- Sanz-Alcaine et al., "Estimation of Semiconductor Power Losses Through Automatic Thermal Modeling", IEEE TPEL 2024, doi 10.1109/TPEL.2024.3412009.
- Repository: github.com/jm-sanz/Automatic-Thermal-Modeling (MIT; `Automatic-Thermal-Modeling-main.zip` MD5 1803fd5f369c7b16c765bfdbe4423f2c).
- Content: 8 thermocouple rises and 6 measured power sources, stored at 100 ms (thermocouples logged at 1 s).
- Used here: the saturation experiment (73 h). The high-side and low-side MOSFETs (Q_h, Q_l) are powered in 2 h blocks with a known power ratio between 0.35 and 5.8.
- Converted to 1 s means in `02_Data/working/sanz_1s/sanz_1s.npz`.

## Coupling is delayed, not flat (`coupling_fast_slow_output.txt`, raw data)
Single-source power step on Q_h (Q_l off), ratio of the Q_l rise to the Q_h rise:

| time after step | 0 s | 5 s | 30 s | 120 s | 600 s | settled |
|---|---|---|---|---|---|---|
| ratio | 0.09 | 0.22 | 0.43 | 0.60 | 0.68 | 0.69 |

- A second step (Q_h off, Q_l on) shows the same pattern: 0.04 at 5 s, rising to 0.41.
- The settled coupling coefficient, from a regression over 35 steps, is 0.59 to 0.66.
- The fast regression is ill-conditioned: the steps are mostly opposite in sign. It is not used.

## Drift estimators under this real coupling (`twin_drift_output.txt`)
- The linear thermal model is identified from the saturation data, as in the source paper.
- Held-out open-loop check (fit 0 to 50 h, simulate 50 to 73 h): RMSE 0.46 K (Q_h) and 0.65 K (Q_l), with rises up to 72 K.
- Q_h is driven with rho x p(t) and Q_l with p(t), where p(t) is the squared mission current of the traction rig. Noise 0.1 K.

| true drift | transfer drift | slope-jump drift |
|---|---|---|
| 0.80 | 0.899 (+12.4 %) | 0.883 to 0.895 (+10 to +12 %) |
| 0.90 | 0.951 (+5.7 %) | 0.939 to 0.949 (+4 to +5 %) |
| 1.10 | 1.046 (−4.9 %) | 1.041 to 1.052 (−4 to −5 %) |
| 1.25 | 1.111 (−11.1 %) | 1.124 to 1.128 (−10 %) |

## Reading
- On a compact board the crosstalk is strong (settled 0.6 to 0.7) and fast (0.2 at 5 s). Both drift estimators are compressed by about one half towards 1.
- The 1 s slope jump is not fast enough to escape the crosstalk on this layout.
- The coupling-free limit exists only for t -> 0: the ratio is 0.09 at the first sample.
- Separate modules on a heat sink are farther apart than two TO devices on one board. How strong the coupling is on the traction rig is not known.
- The audit's fixed-coupling objection is therefore confirmed in size on real hardware. Its flat-coupling form is not: the real coupling grows with time scale.
