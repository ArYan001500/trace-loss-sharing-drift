# E36 Post hoc checks of method v5 (after locked tests L5C and L5F; not locked evidence)

Method v5 and the E29 runners are used unchanged. These checks motivated method v6 (`03_Protocols/METHOD_v6.md`).

## Label swap (`label_swap.py`; outputs `label_swap_{cauer,sink,base}_full.txt`)
- The locked L5 runs are repeated with their locked seeds (Cauer 2032, finite element 2031). Each call of v5 is also evaluated with modules 1 and 2 exchanged in both periods.
  - The interface fault then falls on the reference module, and the true drift is inverted.
- As locked, the L5 results are reproduced exactly.
- With swapped labels:

| run | faults abstained | faults reported beyond 3 % | errors of those (%) | clean abstained |
|---|---|---|---|---|
| Cauer | 0/20 | 20 | −8.6 to −18.6 | 1/40 |
| sink | 0/18 | 2 | +3.4, +3.4 | 1/42 |
| base | 0/18 | 18 | −10.3 to −12.5 | 0/42 |

- Cause: both v5 guards are one-sided, and a fault of module 2 moves both statistics the other way.

## Coupling in the guard windows (`guard_window_coupling.py`; E28 twin)
- Neighbor/own ratio of the slope jumps over the windows of the guards, after a loss step in module 1.
- Maximum over gaps 25/40/55 mm, h 1000/8000 W/m²K, and both placements: [1, 4) 1.2e-5, [1, 5) 6.6e-5, [5, 9) 3.2e-3.
- The same computation on the corrected twin is in E37.

## Assumption of a common loss waveform (`loss_mix_sensitivity.py`; E28 twin, gap 40 mm, h 1000 W/m²K, sensor below the grease, noise-free)
- Each module's loss is a mix of a current-squared and a current-proportional term (mix m1, m2). The loss laws do not change between the periods.
- Different mission phases in the two periods, m up to 0.6: |error| <= 0.37 %.
- Load level of the current period scaled by s (loss laws fixed, so the mix of the two terms changes):
  - s = 0.6: true drift 1.11 to 1.37, errors −3.3 to −6.6 %;
  - s = 0.3: true drift 1.40 to 2.29, errors −9.2 to −14.8 %.
- Reading: the estimand is the change of the loss ratio under a comparable load distribution. A large change of the load level between the periods violates the proportionality assumption.

## Grease resistance of the twin (`fem_twin_grease.m`, `run_fem_grease.m`, `fem_Rg_*.csv`, `fem_grease_compare.py`)
- The E28 twin used R'' = 2.7e-4 m²K/W. The datasheet value of 22.5 K/kW per half module times the half-module footprint (86 mm x 70 mm) gives 1.35e-4 m²K/W.
- Settled coupling below / above the grease, gap 25 mm, h 1000 W/m²K: 0.242 → 0.245 and 0.167 → 0.200. All values stay below c_max = 0.25.
- Interface factor 1.3, gap 40 mm, sensor below the grease, ratio at 4 s: 0.859 → 0.913.
- The full rebuild with the corrected value is E37.

Commands (ROOT = extracted Zenodo record 20395811):
```
python label_swap.py ROOT cauer | sink | base
python guard_window_coupling.py
python loss_mix_sensitivity.py ROOT
matlab -batch "run_fem_grease"
python fem_grease_compare.py
```
