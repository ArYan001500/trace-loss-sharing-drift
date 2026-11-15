# E37 Finite-element twin with the corrected grease resistance

## Why
- The finite-element twin of E28 set the grease area resistance to R'' = 2.7e-4 m²K/W.
- The datasheet gives a case-to-sink resistance of 22.5 K/kW per half module.
  - Multiplied by the half-module footprint (86 mm x 70 mm), this gives R'' = 0.0225 x 0.086 x 0.070 = 1.35e-4 m²K/W, half the E28 value.
  - E28 had in effect applied 22.5 K/kW to the full footprint (86 mm x 140 mm).
  - The comment in `run_all.m` says "over the full module footprint"; the product it computes uses the half-module footprint, as stated here.
- All 24 step responses are rebuilt with this value. Geometry, mesh, sensor points and outputs are unchanged.

## Files
- `fem_twin_grease.m`: E28 `fem_twin.m` with the grease area resistance as the seventh argument.
  - Its header and default (2.7e-4) are those of E28.
  - `run_all.m` always passes 1.35e-4.
  - The file is part of the v6 freeze (`03_Protocols/FREEZE_v6_SHA256.txt`) and was not edited after the freeze.
- `run_all.m`: builds `fem_g{gap}_h{h}_rb{1|1.3}_src{1|2}.csv` for gap 25/40/55 mm, h 1000/8000 W/m²K, interface factor of module 1 1/1.3, and loss step in module 1 or 2.
  - Unlike E28, the gap-40 files also carry the gap in their name.
  - Log: `run_all_output.txt`.
- `fem_field_grease.m`, `fem_field_h1000.csv`, `fem_field_output.txt`: surface field at 4 s, 30 s and 900 s after a loss step in module 1 (gap 40 mm, h 1000 W/m²K), for the figure of the twin.
- `guard_window_coupling.py`, `guard_window_coupling_output.txt`: coupling seen by the guard windows (same computation as E36 on the E28 twin).
- `coupling_bound.py`, `fin_bound.py` and their outputs: copies of E33 run on this twin (plane-source and long-fin screening estimates). `coupling_bound.py` uses the slope law exp(-1/(4 Fo)) of the manuscript; the E33 original used erfc, which is the flux ratio.
- `mesh_check.m`, `mesh_h6.csv`, `mesh_h4.csv`, `mesh_check_output.txt`: maximum element size 6 mm against 4 mm (h 8000, gap 40 mm). Own rise at 900 s 16.894 and 16.918 mK/W (0.14 %); largest difference over 1 to 900 s 0.49 %; settled coupling 0.0224 and 0.0223.
- `fv_crosscheck.py`, `fv_crosscheck_output.txt`: finite-volume model (`src/twin3d.py`, uniform bottom cooling) with the same grease value. Settled coupling 0.174 against 0.193 (h 1000) and 0.016 against 0.022 (h 8000); largest difference 0.020.

Commands (MATLAB R2024b with PDE Toolbox; Python 3):
```
matlab -batch "run_all"
matlab -batch "fem_field_grease(1000, 6e-3, 'fem_field_h1000.csv')"
python guard_window_coupling.py
```

## Results
Settled coupling (neighbor rise / own rise at 900 s), sensor below / above the grease:

| gap (mm) | h = 1000 W/m²K | h = 8000 W/m²K |
|---|---|---|
| 25 | 0.245 / 0.200 | 0.040 / 0.027 |
| 40 | 0.193 / 0.158 | 0.022 / 0.015 |
| 55 | 0.153 / 0.125 | 0.013 / 0.009 |

- The largest value, 0.245, stays below c_max = 0.25 of methods v5 and v6. E28 gave 0.242 at the same point.
- Interface factor 1.3 under module 1, gap 40 mm, own response faulty/clean:
  - sensor below the grease: 0.913 at 4 s, 1.002 at 900 s (E28: 0.859 and 1.003);
  - sensor above the grease, h 1000: 1.104 at 4 s, 1.062 at 900 s (E28: 1.088 and 1.105).
  - The corrected twin therefore shows a smaller interface effect than E28.
- Rise coupling, sensor below the grease, all gaps and cooling levels: at most 3.0e-5 up to 4 s, at most 0.0037 at 10 s, settled 0.013 to 0.245.
- Plane-source estimate: above the twin at 4, 10, 20, 30 and 60 s in all six variants (`coupling_bound_output.txt`).
- Long-fin estimate of the settled coupling: holds in all six variants (`fin_bound_output.txt`).
- Coupling in the guard windows, maximum over gaps, cooling and placements:
  - [1, 4): 2.0e-5; [1, 5): 9.9e-5; [5, 9): 4.7e-3.
  - E28 gave 1.2e-5, 6.6e-5 and 3.2e-3.
