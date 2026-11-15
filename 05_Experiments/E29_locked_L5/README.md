# E29 Locked tests L5C and L5F of method v5

Run once each, as specified in `03_Protocols/METHOD_v5.md` (hashes in `03_Protocols/FREEZE_v5_SHA256.txt`, 41 files).

Development runs made before the freeze:
- `dev_fem_sink_seed11.txt` and `dev_fem_base_seed11.txt` (FEM, seed 11);
- `smoke_cauer_dev7.csv` (Cauer generator, seed 7).

## L5C, Cauer generator (seed 2032, 60 scenarios; `locked_L5C_output.txt`, `locked_L5C_table.csv`)

| criterion | result | required | verdict |
|---|---|---|---|
| D1 clean, abs. error <= 3 % | 1.000 (40/40) | >= 0.90 | PASS |
| D2 interval coverage, clean | 1.000 | >= 0.85 | PASS |
| D3 interface: abstained or <= 3 % | 0.900 (18/20) | >= 0.80 | PASS |
| D4 clean abstention | 0.000 | <= 0.10 | PASS |

Clean error: median 0.41 %, maximum 2.13 %.

## L5F, finite-element twin (seed 2031, 60 scenarios per sensor placement)

| placement | D1 | D2 | D3 | D4 | clean median / max |
|---|---|---|---|---|---|
| heat-sink surface, below the grease (`locked_L5F_sink_output.txt`) | 1.000 (41/41) | 1.000 | 1.000 (18/18) | 0.024 | 0.53 / 2.82 % |
| baseplate bottom, above the grease (`locked_L5F_base_output.txt`) | 1.000 (39/39) | 1.000 | 1.000 (18/18) | 0.071 | 0.20 / 0.89 % |

- The draws include module gaps of 25, 40 and 55 mm, cooling h of 1000 and 8000 W/m²K, and the thermal coupling between modules.
- Clean abstentions on the base placement: 3 of 42, all from the shape guard (z −1.97 to −2.74). Their errors would have been below 0.4 %.

## Verdict
Method v5 passes all four criteria in all three locked runs.

## Limits
- Both generators are simulations. The finite-element geometry is estimated from figures.
- No real recording with a known loss-ratio change has been tested.
