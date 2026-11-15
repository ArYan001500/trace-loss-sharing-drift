"""Figures for the Results section, drawn at final printed size (no scaling in LaTeX).
Column width 86.4 mm and page width 177 mm (the template text block), Times New Roman 8 pt, vector PDF with embedded TrueType fonts.
python make_results_figs.py PROJECT_ROOT OUT_DIR"""
import sys, os, re, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT, OUT = sys.argv[1], sys.argv[2]
EXP = os.path.join(ROOT, "05_Experiments")
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(EXP, "E05_mission_derivatives"))
MM = 1 / 25.4
COL, PAGE = 86.4 * MM, 177 * MM
TWIN_H, RIG_H = float(os.environ.get("TWIN_H", 46)), float(os.environ.get("RIG_H", 42))   # figure heights in mm
plt.rcParams.update({
    "font.family": "Times New Roman", "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.0,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False})
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
          blue="#0072B2", red="#D55E00", purple="#CC79A7", grey="#888888")


def panel_label(ax, s, x=-0.16, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches=None)
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- Fig. twin: field + coupling
def fig_twin():
    fem = os.path.join(EXP, "E37_fem_twin_grease")
    fig = plt.figure(figsize=(PAGE, TWIN_H * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], left=0.075, right=0.99, bottom=10.2 / TWIN_H, top=1 - 4.8 / TWIN_H, wspace=0.72)
    F = np.loadtxt(os.path.join(fem, "fem_field_h1000.csv"), delimiter=",")
    x, y = F[:, 0] * 1e3, F[:, 1] * 1e3
    nx, ny = len(np.unique(x)), len(np.unique(y))
    X, Y = x.reshape((ny, nx), order="F"), y.reshape((ny, nx), order="F")
    mw, ml, y0, x1, x2 = 86, 140, 40, 30, 156
    for j, (col, t, lab) in enumerate(((2, 4, "(a)"), (4, 900, "(b)"))):
        ax = fig.add_subplot(gs[0, j])
        Z = F[:, col].reshape((ny, nx), order="F") * 1e3          # rise in K per kW of loss in module 1
        im = ax.pcolormesh(X, Y, Z, cmap="inferno", shading="auto", rasterized=True)
        for xx in (x1, x2):
            ax.add_patch(Rectangle((xx, y0), mw, ml, fill=False, ec="white", lw=0.7, ls="--"))
        for xx in (x1 + mw / 2, x2 + mw / 2):
            for f in (0.69, 0.64, 0.33, 0.28):
                ax.plot(xx, y0 + f * ml, "o", ms=2.2, mfc=OI["sky"], mec="white", mew=0.4)
        ax.text(x1 + mw / 2, y0 + ml + 8, "Module 1", color="white", ha="center")
        ax.text(x2 + mw / 2, y0 + ml + 8, "Module 2", color="white", ha="center")
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        ax.set_xlim(0, X.max()); ax.set_ylim(0, Y.max())
        ax.spines["top"].set_visible(True); ax.spines["right"].set_visible(True)
        ax.set_title(f"{t} s after the loss step", pad=3)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.set_label("Rise (K/kW)", labelpad=2); cb.outline.set_linewidth(0.5)
        panel_label(ax, lab, x=-0.3)
    ax = fig.add_subplot(gs[0, 2])
    t = np.arange(1, 901)
    gc = {25: OI["red"], 40: OI["orange"], 55: OI["blue"]}
    styles = {(g, h): ("-" if h == 1000 else "--", gc[g]) for g in (25, 40, 55) for h in (1000, 8000)}
    for (gp, h), (ls, c) in styles.items():
        fn = f"fem_g{gp}_h{h}_rb1_src1.csv"
        S = np.loadtxt(os.path.join(fem, fn), delimiter=",")[:, :8]
        with np.errstate(divide="ignore", invalid="ignore"):
            cpl = S[1:, 4:].mean(1) / S[1:, :4].mean(1)
        ax.plot(t, np.clip(cpl, 1e-5, None), ls=ls, color=c, lw=1.0, label=(f"Gap {gp} mm" if h == 1000 else None))
    # measured board: neighbor rise / own rise after a loss step in one transistor
    d = np.load(os.path.join(ROOT, "02_Data", "working", "sanz_1s", "sanz_1s.npz"))["Sat"]
    s = 7201; tt = np.arange(1, 901)
    cb_ = (d[1, s + tt] - d[1, s - 1]) / (d[0, s + tt] - d[0, s - 1])
    ax.plot(tt, cb_, color=OI["black"], lw=1.2, label="Board")
    a = 167 / (2700 * 900)
    from scipy.special import erfc
    eta = 0.083 / (2 * np.sqrt(a * t))
    ax.plot(t, np.exp(-eta ** 2) - np.sqrt(np.pi) * eta * erfc(eta), ":", color=OI["grey"], lw=1.2, label="Plane source")
    ax.axvspan(1, 4, color=OI["grey"], alpha=0.15, lw=0)
    ax.axhline(0.25, color=OI["purple"], lw=0.9, ls=(0, (4, 2)), zorder=1)
    ax.text(1.12, 0.3, "$c$" + "max = 0.25", color=OI["purple"], ha="left", va="bottom", fontsize=8)
    ax.text(2.0, 1.5e-3, "early\nwindow", ha="center", va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-5, 1.5); ax.set_xlim(1, 900)
    ax.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]); ax.set_yticklabels(["0.00001", "0.0001", "0.001", "0.01", "0.1", "1"])
    ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
    ax.minorticks_off()
    ax.set_xlabel("Time after the loss step (s)"); ax.set_ylabel("Neighbor rise / own rise")
    ax.legend(loc="lower right", handlelength=1.6, handletextpad=0.4, borderaxespad=0.1, labelspacing=0.15, frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    panel_label(ax, "(c)", x=-0.3)
    save(fig, "fig_twin.pdf")


# ---------------------------------------------------------------- Fig. locked tests
def parse_fem_locked(path):
    rows = []
    for l in open(path, encoding="utf-8"):
        m = re.match(r"\s*\d+ gap (\d+) h (\d+) fault (\d) truth ([\d.]+) est ([\d.]+) \(([-+\d.]+) %\)(.*)", l)
        if m:
            rows.append(dict(fault=int(m.group(3)), err=float(m.group(6)), abst="ABSTAIN" in m.group(7)))
    return pd.DataFrame(rows)


def fig_locked():
    L = os.path.join(EXP, "E38_locked_L6")
    def tab(f):
        c = pd.read_csv(os.path.join(L, f))
        return pd.DataFrame(dict(fault=c.fault.astype(int), err=c.err_pct, abst=c.abstain))
    runs = [("Datasheet\ngenerator", tab("locked_L6C_table.csv")),
            ("Twin, sensor\nbelow grease", tab("locked_L6F_sink_table.csv")),
            ("Twin, sensor\nabove grease", tab("locked_L6F_base_table.csv"))]
    fig, ax = plt.subplots(figsize=(COL, 58 * MM))
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.19, top=0.86)
    rng = np.random.default_rng(0)
    ax.axhspan(-3, 3, color=OI["grey"], alpha=0.15, lw=0)
    for i, (name, r) in enumerate(runs):
        for fault, (mk, col, lab) in {0: ("o", OI["blue"], "Clean cases"), 1: ("^", OI["red"], "Interface faults")}.items():
            q = r[(r.fault == fault) & (~r.abst)]
            off = -0.17 if fault == 0 else 0.17
            ax.scatter(i + off + rng.uniform(-0.09, 0.09, len(q)), q.err, s=10, marker=mk, color=col, lw=0,
                       label=lab if i == 0 else None, zorder=3)
            tot = (r.fault == fault).sum(); ab = (r.abst & (r.fault == fault)).sum()
            ax.text(i + off, 14.0, f"{ab}/{tot}", ha="center", color=col)
        ax.text(i, 16.4, "abstained", ha="center", color=OI["black"]) if i == 0 else None
    ax.set_xticks(range(3)); ax.set_xticklabels([n for n, _ in runs])
    ax.set_ylabel("Drift error of reported cases (%)"); ax.set_ylim(-15, 19); ax.set_xlim(-0.5, 2.5)
    ax.axhline(0, color="black", lw=0.4)
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.2), columnspacing=1.0, handletextpad=0.2)
    save(fig, "fig_locked.pdf")


# ---------------------------------------------------------------- Fig. comparison
def parse_cmp(path):
    rows = []
    for l in open(path, encoding="utf-8"):
        m = re.match(r"\s*\d+ fault (\d) side (\d) truth ([\d.]+) \| (.*?)(\| v6 abstains)?\s*$", l)
        if m:
            vals = dict(re.findall(r"(static|wavelet|twotemp|slope|v6|calib)\s+([-+\d.]+)", m.group(4)))
            rows.append(dict(fault=int(m.group(1)), side=int(m.group(2)), abst=m.group(5) is not None,
                             **{k: float(v) for k, v in vals.items()}))
    return pd.DataFrame(rows)


def fig_compare():
    E = os.path.join(EXP, "E39_v6_support")
    names = [("static", "Mean-rise ratio"), ("wavelet", "Wavelet ratio"), ("twotemp", "Two-temperature fit"),
             ("slope", "Slope jump"), ("calib", "Calibrated inversion*"), ("v6", "TRACE")]
    fig, axes = plt.subplots(1, 2, figsize=(COL, 65 * MM), sharey=True)
    fig.subplots_adjust(left=0.335, right=0.985, bottom=0.16, top=0.775, wspace=0.1)
    for ax, (fn, title) in zip(axes, (("cmp_cauer_output.txt", "Generator"), ("cmp_fem_output.txt", "Twin"))):
        r = parse_cmp(os.path.join(E, fn))
        ax.axhspan(-0.5, 0.5, color="#FFF0C2", lw=0, zorder=0)          # highlight the TRACE row
        ax.axvspan(-3, 3, color=OI["grey"], alpha=0.15, lw=0)
        for j, (k, lab) in enumerate(names):
            yv = len(names) - 1 - j
            sel = (r.fault == 0) & (~r.abst if k == "v6" else True)        # TRACE: reported clean cases only
            cl = r[sel][k].clip(-20, 22)
            ax.scatter(cl, yv + np.random.default_rng(j).uniform(-0.18, 0.18, len(cl)), s=5, color=OI["blue"], alpha=0.75, lw=0, zorder=3)
            for sd, mk in ((1, "^"), (2, "v")):
                g = r[(r.fault == 1) & (r.side == sd)]
                im = g[k].median()
                hollow = k == "v6" and bool(g.abst.all())                  # hollow only if TRACE abstained in every case
                ax.plot(np.clip(im, -20, 22), yv, marker=mk, ms=5, mfc="white" if hollow else OI["red"], mec=OI["red"], mew=1.0, zorder=4)
        ax.set_title(title, pad=3); ax.set_xlim(-22, 24); ax.set_xticks([-20, -10, 0, 10, 20]); ax.set_xlabel("Drift error (%)")
        ax.set_ylim(-0.5, len(names) - 0.5)
        ax.axvline(0, color="black", lw=0.4)
    axes[0].set_yticks(range(len(names))); axes[0].set_yticklabels([lab for _, lab in names][::-1])
    axes[0].get_yticklabels()[0].set_fontweight("bold")
    h1, = axes[0].plot([], [], "o", ms=3, color=OI["blue"])
    h2, = axes[0].plot([], [], "^", ms=5, color=OI["red"])
    h4, = axes[0].plot([], [], "v", ms=5, color=OI["red"])
    h3, = axes[0].plot([], [], "^", ms=5, mfc="white", mec=OI["red"])
    fig.legend([h1, h2, h4, h3], ["Clean cases", "Fault, module 1", "Fault, module 2", "All abstained"], loc="upper center", ncol=2,
               bbox_to_anchor=(0.62, 1.0), columnspacing=0.8, handletextpad=0.1, handlelength=1.0, labelspacing=0.2)
    save(fig, "fig_compare.pdf")


# ---------------------------------------------------------------- Fig. rig case study
def fig_rig():
    from load_mission import load
    import tdrift, method_v6
    from wcoh import coherence
    root = os.path.join(ROOT, "02_Data", "raw", "zenodo_20395811", "extracted")
    df = load(root, 0); n = len(df); h = n // 2
    PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}; CIRC = {"U": "A", "V": "A", "W": "C"}
    fig = plt.figure(figsize=(PAGE, RIG_H * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 0.8, 1.05], left=0.068, right=0.935, bottom=11.4 / RIG_H, top=1 - 6.0 / RIG_H, wspace=0.9)
    ax = fig.add_subplot(gs[0, 0]); t = np.arange(n)
    I = df["Ia_ef"].shift(-11).to_numpy()
    ax.plot(t, df["U1T4"], color=OI["red"], lw=0.9, label="Module 1")
    ax.plot(t, df["U2T4"], color=OI["blue"], lw=0.9, label="Module 2")
    Tw = 0.5 * (df["Tw_A-In"] + df["Tw_A-Out"])
    ax.plot(t, Tw, color=OI["green"], lw=0.9, label="Coolant")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Temperature (°C)"); ax.set_xlim(0, n)
    ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
    ax2.fill_between(t, 0, I, color=OI["grey"], alpha=0.25, lw=0, label="Current")
    ax2.set_ylabel("Phase current (A)"); ax2.set_ylim(0, 2600)
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    ax.axvline(h, color="black", lw=0.5, ls="--"); ax.text(h + 6, 26.5, "halves")
    l1, b1 = ax.get_legend_handles_labels(); l2, b2 = ax2.get_legend_handles_labels()
    ax.legend(l1 + l2, b1 + b2, loc="upper left", handlelength=1.2, labelspacing=0.15, ncol=2, columnspacing=0.6, borderaxespad=0.1)
    ax.set_ylim(24, 60)
    panel_label(ax, "(a)", x=-0.13, y=1.05)
    ax = fig.add_subplot(gs[0, 1])
    for i, ph in enumerate("UVW"):
        Tc = (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()
        pairs = [(f"T{k}", df[f"{ph}1T{k}"].to_numpy(), df[f"{ph}2T{k}"].to_numpy()) for k in "1234" if f"{ph}1T{k}" in df and f"{ph}2T{k}" in df]
        cut = lambda s: [(nm, a, b) for nm, a, b in [(nm, a[s], b[s]) for nm, a, b in pairs]]
        Ip = df[PH[ph]].shift(-11).ffill().bfill().to_numpy()
        r = method_v6.estimate_drift((Ip[:h], Tc[:h], cut(slice(0, h))), (Ip[h:], Tc[h:], cut(slice(h, n))))
        td = tdrift.estimate_drift(Tc[:h], cut(slice(0, h)), Tc[h:], cut(slice(h, n)))
        pos = [np.exp(p["log_d"]) for p in td["positions"].values()]
        ax.plot([i] * len(pos), [100 * (p - 1) for p in pos], "o", ms=3, mfc="white", mec=OI["blue"], mew=0.8,
                label="Positions" if i == 0 else None)
        lo, hi = r["ci95"]
        ax.errorbar(i + 0.22, 100 * (r["drift"] - 1), yerr=[[100 * (r["drift"] - lo)], [100 * (hi - r["drift"])]],
                    fmt="s", ms=3.5, color=OI["red"], capsize=2, lw=0.9, label="TRACE" if i == 0 else None)
    ax.axhline(0, color="black", lw=0.4)
    ax.set_xticks(range(3)); ax.set_xticklabels(["U", "V", "W"]); ax.set_xlim(-0.5, 2.6); ax.set_xlabel("Phase")
    ax.set_ylabel("Drift (%)"); ax.set_ylim(-2.0, 3.0)
    ax.legend(loc="upper center", ncol=2, handletextpad=0.1, columnspacing=0.6, borderaxespad=0.1)
    panel_label(ax, "(b)", x=-0.38, y=1.05)
    ax = fig.add_subplot(gs[0, 2])
    base = os.path.join(ROOT, "02_Data", "raw", "zenodo_20395811", "extracted", "Parallel modules under different cooling circuits MAG_1_1")
    def run(name):
        x = pd.read_excel(os.path.join(base, name, f"MAG1_1_{name}.xlsx"))
        x.columns = ["time"] + [str(c).split(" - ")[-1].strip() for c in x.columns[1:]]
        tt = (pd.to_datetime(x.time) - pd.to_datetime(x.time[0])).dt.total_seconds().to_numpy()
        return x, tt
    periods = np.array([40.0, 60, 90, 130, 200, 300]); res = {}
    for name in ("L0_0", "L0_50"):
        x, tt = run(name); tg = np.arange(0, tt[-1], 10.0)
        Tc = np.interp(tg, tt, x["Tw_C-In"].to_numpy(float)); res[name] = {}
        for k in "1234":
            a = np.interp(tg, tt, x[f"W1T{k}"].to_numpy(float)) - Tc; b = np.interp(tg, tt, x[f"W2T{k}"].to_numpy(float)) - Tc
            R2, K, P = coherence(a, b, dt=10.0, periods=periods)
            m = len(tg); e = slice(int(0.15 * m), int(0.85 * m))
            res[name][k] = (np.median(np.abs(K[:, e]), axis=1), np.degrees(np.median(np.angle(K[:, e]), axis=1)))
    ax3 = ax.twinx(); ax3.spines["right"].set_visible(True)
    cols = [OI["blue"], OI["sky"], OI["green"], OI["orange"]]
    for c, k in zip(cols, "1234"):
        ax.plot(periods, 100 * (res["L0_50"][k][0] / res["L0_0"][k][0] - 1), "o-", ms=2.5, color=c, lw=0.9, label=f"T{k}")
        ax3.plot(periods, res["L0_50"][k][1] - res["L0_0"][k][1], "x:", ms=3, color=c, lw=0.7)
    ax.axhline(0, color="black", lw=0.4); ax.set_xscale("log")
    ax.set_xticks([40, 60, 100, 200, 300]); ax.set_xticklabels(["40", "60", "100", "200", "300"])
    ax.set_xlabel("Period (s)"); ax.set_ylabel("Ratio change (%)"); ax.set_ylim(-10, 4)
    ax3.set_ylabel("Phase change (deg)"); ax3.set_ylim(-5, 9); ax3.set_yticks(range(-4, 9, 2)); ax.set_yticks(range(-10, 5, 2))
    ax.legend(loc="upper right", ncol=4, columnspacing=0.5, handlelength=1.0, handletextpad=0.2, borderaxespad=0.1)
    panel_label(ax, "(c)", x=-0.3, y=1.05)
    save(fig, "fig_rig.pdf")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in (fig_twin, fig_locked, fig_compare, fig_rig):
        f()
