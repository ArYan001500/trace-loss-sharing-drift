# E38 Locked tests L6C and L6F of method v6

Protocol: `03_Protocols/METHOD_v6.md`, frozen with `03_Protocols/FREEZE_v6_SHA256.txt` before the locked runs.
The hashes were verified (`sha256sum -c`) immediately before the runs.

## Files
- `run_v6_cauer.py`, `run_v6_fem.py`, `common.py`: runners and the summary of criteria D1 to D4.
- `dev_*`: development runs (Cauer seed 7, finite-element seed 11, 60 scenarios each), made before the freeze.
- `locked_L6C_*`, `locked_L6F_sink_*`, `locked_L6F_base_*`: locked runs (seeds 2041 and 2042, 100 scenarios each), each run once.

Commands (ROOT = extracted Zenodo record 20395811):
```
python run_v6_cauer.py ROOT 2041 100 locked_L6C_table.csv 3.0
python run_v6_fem.py ROOT sink 2042 100 locked_L6F_sink_table.csv
python run_v6_fem.py ROOT base 2042 100 locked_L6F_base_table.csv
```

## Locked results (appended after the runs; code and protocol unchanged)

| run | clean n / reported | interface n (module 1 / 2) | D1 | D2 | D3 | D4 | result |
|---|---|---|---|---|---|---|---|
| L6C Cauer | 81 / 81 | 19 (13 / 6) | 1.000 | 1.000 | 0.895 | 0.000 | PASS |
| L6F sink | 68 / 65 | 32 (17 / 15) | 0.969 | 1.000 | 1.000 | 0.044 | PASS |
| L6F base | 68 / 65 | 32 (17 / 15) | 1.000 | 1.000 | 0.812 | 0.044 | PASS |

Thresholds: D1 >= 0.90, D2 >= 0.85, D3 >= 0.80, D4 <= 0.10. Method v6 passes all three runs.

D3 counts per faulted module (abstained / reported within 3 % / reported beyond 3 %):
- L6C: module 1, 13 / 0 / 0; module 2, 4 / 0 / 2 (errors −12.9 % and −9.6 %).
- L6F sink: module 1, 17 / 0 / 0; module 2, 15 / 0 / 0.
- L6F base: module 1, 13 / 0 / 4; module 2, 13 / 0 / 2 (errors 8.5 % to 11.6 % in magnitude).
  - In these six cases the fault moved the fast and the wideband drift together (|z_x| <= 2.14).
  - The sensor above the grease sees the interface change partly as a gain change, which the guards cannot separate from a loss change.

Clean reported cases:
- Median |error|: 0.52 % (Cauer), 0.47 % (sink), 0.24 % (base).
- Maximum |error|: 2.39 % (Cauer), 4.86 % (sink), 2.58 % (base).
- The two sink cases beyond 3 % (−4.86 % and +3.05 %) are at a module gap of 25 or 40 mm and h = 1000 W/m²K, where settled coupling compresses the wideband drift.
  - Both lie inside the coupling-extended interval.

Clean abstentions: 0 (Cauer), 3 (sink), 3 (base). Codes: `band_disagreement` (4) and `shape_change` (2).

## Relation to earlier locked tests
- L5C and L5F (method v5, E29) passed with faults on module 1 only.
- A post hoc label swap (E36) showed that v5 missed every module-2 fault.
- L6 draws the faulted module at random and uses the finite-element twin with the corrected grease resistance (E37).
- The L5F and L6F base results are therefore not comparable: the grease resistance differs by a factor of two.
