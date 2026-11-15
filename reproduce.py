"""Downloads the public data and reproduces the results of the paper. Outputs go to results/.
python reproduce.py            all steps (about one hour)
python reproduce.py --locked   download, hash check and the three locked tests only"""
import hashlib, os, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
R = os.path.join("02_Data", "raw", "zenodo_20395811", "extracted")
NPZ = os.path.join("02_Data", "working", "sanz_1s", "sanz_1s.npz")
E = "05_Experiments"
PY = sys.executable


def run(args, out, cwd=ROOT):
    print(">", " ".join(args[1:]), flush=True)
    with open(os.path.join(RES, out), "w", encoding="utf-8") as f:
        subprocess.run(args, cwd=cwd, stdout=f, stderr=subprocess.DEVNULL, check=True)


def freeze_check():
    bad = []
    for line in open(os.path.join(ROOT, "03_Protocols", "FREEZE_v6_SHA256.txt"), encoding="utf-8"):
        digest, path = line.split(maxsplit=1)
        path = path.strip().lstrip("*")
        with open(os.path.join(ROOT, path), "rb") as f:
            if hashlib.sha256(f.read()).hexdigest() != digest:
                bad.append(path)
    if bad:
        sys.exit("files changed since the freeze: " + ", ".join(bad))
    print("all files of the method v6 freeze match their SHA-256 hashes")


def main(locked_only):
    os.makedirs(RES, exist_ok=True)
    subprocess.run([PY, os.path.join("data", "download_data.py")], cwd=ROOT, check=True)
    freeze_check()

    # Locked tests L6 (Table 2, Fig. 2), compared row by row with the stored tables
    run([PY, f"{E}/E38_locked_L6/run_v6_cauer.py", R, "2041", "100", os.path.join(RES, "locked_L6C_table.csv"), "3.0"], "locked_L6C_output.txt")
    for p in ("sink", "base"):
        run([PY, f"{E}/E38_locked_L6/run_v6_fem.py", R, p, "2042", "100", os.path.join(RES, f"locked_L6F_{p}_table.csv")], f"locked_L6F_{p}_output.txt")
    for t in ("L6C", "L6F_sink", "L6F_base"):
        subprocess.run([PY, "verify_tables.py", os.path.join(RES, f"locked_{t}_table.csv"), f"{E}/E38_locked_L6/locked_{t}_table.csv"], cwd=ROOT, check=True)
    if locked_only:
        return

    # Development runs, ablation, comparison, worked example, load-mix probe (Table 3, Fig. 3, Sections III-E, IV-A)
    run([PY, f"{E}/E38_locked_L6/run_v6_cauer.py", R, "7", "60", os.path.join(RES, "dev_cauer_seed7.csv")], "dev_cauer_seed7.txt")
    for g, s in (("cauer", "7"), ("fem_sink", "11"), ("fem_base", "11")):
        run([PY, f"{E}/E39_v6_support/ablation_v6.py", R, g, s, "60"], f"abl_{g}.txt")
    for g in ("cauer", "fem"):
        run([PY, f"{E}/E39_v6_support/compare_baselines_v6.py", R, g, "101", "40"], f"cmp_{g}_output.txt")
    run([PY, f"{E}/E39_v6_support/worked_example_v6.py", R], "worked_example_v6_output.txt")
    run([PY, f"{E}/E39_v6_support/loss_mix_sensitivity.py", R], "loss_mix_sensitivity_output.txt")

    # Rig case study and board negative control (Sections VI-D, VI-E)
    run([PY, f"{E}/E39_v6_support/rig_v6_halves.py", R], "rig_v6_halves_output.txt")
    run([PY, f"{E}/E39_v6_support/rig_sensitivity.py", R], "rig_sensitivity_output.txt")
    run([PY, f"{E}/E23_transfer_drift/rig_cooling_runs.py", os.path.join(R, "Parallel modules under different cooling circuits MAG_1_1")], "rig_cooling_runs_output.txt")
    run([PY, f"{E}/E39_v6_support/board_control.py", NPZ, R, "41", "40"], "board_control_output.txt")
    run([PY, f"{E}/E39_v6_support/board_control_qonly.py", NPZ, R, "41", "40"], "board_control_qonly_output.txt")
    run([PY, f"{E}/E25_real_coupling/coupling_fast_slow.py", NPZ], "board_coupling_output.txt")

    # Coupling of the twin (Section VI-D)
    twin = os.path.join(ROOT, E, "E37_fem_twin_grease")
    run([PY, "guard_window_coupling.py"], "twin_guard_windows.txt", cwd=twin)
    run([PY, "fin_bound.py", "."], "twin_fin_bound.txt", cwd=twin)
    run([PY, "coupling_bound.py", ".", os.path.join(ROOT, NPZ)], "twin_plane_source.txt", cwd=twin)

    # Post hoc checks of method v5 (Table 4) and of the coupling allowance (Section VI-A)
    for r in ("cauer", "sink", "base"):
        run([PY, f"{E}/E36_posthoc_review/label_swap.py", R, r], f"label_swap_{r}.txt")
        run([PY, "-W", "ignore", f"{E}/E40_interval_audit/interval_audit.py", R, r], f"interval_audit_{r}.txt")

    # Figures (Figs. 2 and 3 are drawn from the stored tables in E38 and E39, which the steps above reproduce)
    fig = os.path.join(RES, "figures"); os.makedirs(fig, exist_ok=True)
    subprocess.run([PY, "07_Manuscript/figures_src/make_results_figs.py", ".", fig], cwd=ROOT, check=True, stderr=subprocess.DEVNULL)
    subprocess.run([PY, "07_Manuscript/figures_src/make_overview_fig_3d.py", fig], cwd=ROOT, check=True)
    print("done; outputs in", RES)


if __name__ == "__main__":
    main("--locked" in sys.argv)
