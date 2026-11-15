# TRACE: temperature-based monitoring of loss-sharing drift in parallel power modules

Code and data for the paper "Temperature-Based Monitoring of Loss-Sharing Drift in Parallel Power Modules With Consistency-Based Abstention" by Aryan Ashoori.

TRACE estimates how the loss sharing of two parallel power modules changes against a healthy baseline, using module and coolant temperatures. The phase current is used only to locate load steps. Two response-consistency checks withhold the estimate when the thermal responses do not fit a pure change of loss sharing.

## Quick start

```
pip install -r requirements.txt
python reproduce.py
```

`reproduce.py` downloads the public data, checks the frozen files of the method against their SHA-256 hashes, reruns the locked tests and compares every row with the stored tables, and then reruns the supporting analyses and the figures. A full run takes about one hour. `python reproduce.py --locked` stops after the three locked tests. Outputs are written to `results/`.

## Data

`data/download_data.py` fetches two public datasets. It is called by `reproduce.py`.

- Traction-rig recordings, Zenodo record 20395811 (CC BY 4.0). The script downloads two ZIP archives, checks their MD5 sums, and extracts the mission recording and the two cooling-circuit runs.
- Buck-converter thermal dataset of Sanz-Alcaine et al., github.com/jm-sanz/Automatic-Thermal-Modeling (MIT license), pinned to a fixed commit. The script extracts `Dataset.rar`, checks the four CSV files, and builds 1 s means with `data/make_sanz_1s.py`. Extracting the RAR file needs bsdtar (the `tar` of Windows 10 and later, or `libarchive-tools` on Linux), 7z or unrar.

## Layout

- `src/`: the estimators. `method_v6.py` is TRACE as evaluated in the paper. It uses `tdrift.py` for the transfer-referenced drift and `method_v3.py` with `fingerprint_est.py` for the load steps and the guards. `cauer_ds.py` is the datasheet-based generator, `twotemp.py` and `wcoh.py` are comparison methods, and `twin3d.py` is a finite-volume model used as a cross-check. `method_v5.py` is the previous version, kept for the post hoc checks in the paper.
- `03_Protocols/`: the frozen protocols of versions 5 and 6 and their SHA-256 hash lists. The hashes were recorded locally before each locked test and carry no external time stamp.
- `05_Experiments/`: one folder per study, each with a README and a hash list.

| folder | content | paper |
|---|---|---|
| E37_fem_twin_grease | finite-element twin (MATLAB PDE Toolbox), stored step responses, coupling checks | Sections V-A and VI-D, Fig. 4 |
| E38_locked_L6 | runners, development runs, and the locked tests L6C and L6F with per-case tables | Table 2, Fig. 2 |
| E39_v6_support | ablation, comparison methods, rig, board, worked example, load-mix probe | Table 3, Fig. 3, Sections III-E, IV-A, VI |
| E40_interval_audit | coverage with and without the coupling allowance on the locked scenarios | Section VI-A |
| E36_posthoc_review | label swap and grease check of version 5 | Table 4, Section VII-B |
| E28_fem_twin, E29_locked_L5 | earlier twin and the version 5 runners, inputs of E36 | Table 4 |
| E23_transfer_drift | exploratory comparison of the two cooling runs, logged at 10 s | Section VI-E |
| E25_real_coupling | coupling measured on the board | Section VI-D |
| E32_matlab_crosscheck | an independent MATLAB implementation of the wideband estimator | Section VI-E |
| E05, E16 | mission loader and generator parameters | Section V |

- `07_Manuscript/figures_src/`: the scripts behind every figure.
- `verify_tables.py`: compares a recomputed locked-test table with the stored one, row by row, and recomputes the criteria D1 to D4.

Each locked-test table in `E38_locked_L6` lists, for every scenario, the drawn conditions, the true drift, the estimate, the interval, the guard statistics, the abstention decision and its reason code. The seeds are given in `03_Protocols/METHOD_v6.md`.

## Requirements

Python 3 with the packages in `requirements.txt` (the versions used for the paper). The finite-element step responses are stored, so MATLAB is not needed to reproduce the results. Recomputing them needs MATLAB R2024b with the PDE Toolbox (`05_Experiments/E37_fem_twin_grease/run_all.m`). The MATLAB cross-check in `E32` also needs MATLAB and is not part of `reproduce.py`.

## Contact

Aryan Ashoori, aryanashoori0@gmail.com
