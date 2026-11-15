"""Fig. 1 (candidate): (a) isometric 3-D view of two parallel modules on a liquid-cooled plate with the
heat of module 1 spreading, (b) TRACE procedure as numbered steps.
Drawn at the template text width (177 mm), Times New Roman 8 pt, vector PDF.
python make_overview_fig_3d.py OUT_DIR"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch, FancyArrowPatch, Circle, Ellipse
from matplotlib.colors import LinearSegmentedColormap

OUT = sys.argv[1]
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix", "font.size": 8, "pdf.fonttype": 42})
C = dict(blue="#0072B2", red="#D55E00", green="#009E73", grey="#777777", sky="#56B4E9",
         orange="#E69F00", dark="#3A3A3A", copper="#C77B3A", purple="#9B3F8A")
C30, S30 = np.cos(np.pi / 6), np.sin(np.pi / 6)


def P(x, y, z):
    """Isometric projection: x to the right and up, y to the left and up, z up."""
    return np.array([(x - y) * C30, (x + y) * S30 + z])


def face(ax, pts, fc, ec="black", lw=0.5, z=1, alpha=1.0):
    ax.add_patch(Polygon([P(*p) for p in pts], closed=True, fc=fc, ec=ec, lw=lw, zorder=z, alpha=alpha))


def box(ax, x0, y0, z0, dx, dy, dz, top, left, right, z=1, lw=0.5):
    x1, y1, z1 = x0 + dx, y0 + dy, z0 + dz
    face(ax, [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)], right, lw=lw, z=z)
    face(ax, [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)], left, lw=lw, z=z)
    face(ax, [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)], top, lw=lw, z=z)


fig = plt.figure(figsize=(177 * MM, 62 * MM))

# ============================================================ (a) 3-D setup
ax = fig.add_axes([0.0, 0.0, 0.43, 1.0]); ax.set_aspect("equal"); ax.axis("off")
PX, PY, PZ = 150.0, 62.0, 9.0
box(ax, 0, 0, 0, PX, PY, PZ, top="#DCE6EE", left="#AFC3D3", right="#C3D3E0", z=1)
# heat of module 1 on the plate surface, spreading towards module 2
heat = LinearSegmentedColormap.from_list("heat", ["#FFF1D0", "#F9C46B", "#EE7D2E", "#C8331B"])
cx, cy = 34.0, 31.0
th = np.linspace(0, 2 * np.pi, 160)
for k, r in enumerate(np.linspace(1.0, 0.3, 8)):
    xs = np.clip(cx + 64 * r * np.cos(th), 0.4, PX - 0.4); ys = np.clip(cy + 36 * r * np.sin(th), 0.4, PY - 0.4)
    ax.add_patch(Polygon([P(x, y, PZ) for x, y in zip(xs, ys)], closed=True, fc=heat(k / 7), ec="none",
                         alpha=0.75, zorder=2))
# coolant channel openings and flow
for yy in (15, 31, 47):
    ax.add_patch(Ellipse(P(0, yy, PZ / 2), 5.5, 3.6, angle=-30, fc="#1D4E7A", ec="black", lw=0.4, zorder=3))
ax.add_patch(FancyArrowPatch(P(-34, 31, PZ / 2), P(-5, 31, PZ / 2), arrowstyle="-|>", mutation_scale=9, lw=1.5,
                             color=C["blue"], zorder=4))
q = P(-34, 31, PZ / 2); ax.text(q[0] + 1, q[1] - 2.5, "Coolant", ha="left", va="top", color=C["blue"])
# thermocouples at matched positions (front row shown)
for mx in (12.0, 88.0):
    for k in range(4):
        ax.add_patch(Circle(P(mx + 7 + 10 * k, 2.6, PZ), 1.5, fc=C["orange"], ec="black", lw=0.4, zorder=5))
# modules
mods = [(12.0, "Module 1", C["red"], r"loss $\rho\,p(t)$"), (88.0, "Module 2", C["blue"], r"loss $p(t)$")]
for mx, name, col, loss in mods[::-1]:
    box(ax, mx, 6, PZ, 46, 50, 3.0, top="#9C9C9C", left="#6F6F6F", right="#828282", z=7)
    box(ax, mx + 3, 9, PZ + 3, 40, 44, 12, top="#F5F5F5", left="#BDBDBD", right="#D8D8D8", z=8)
    for yy in (14, 27, 40):
        box(ax, mx + 15, yy, PZ + 15, 10, 7, 3.5, top="#E8B07A", left=C["copper"], right="#D88E4C", z=9)
    c = P(mx + 4, 56, PZ + 15) if mx < 50 else P(PX, 22, PZ + 16)
    ha = "center" if mx < 50 else "left"
    cx_, cy_ = (c[0] - 6, c[1] + 15) if mx < 50 else (c[0] + 1, c[1] + 2)
    ax.text(cx_, cy_, name, ha=ha, va="bottom", color=col, fontweight="bold", zorder=12)
    ax.text(cx_, cy_ - 0.8, loss, ha=ha, va="top", color=col, zorder=12)
# phase current: one label, a junction, and two straight branches to the terminals
J = P(54, 31, PZ + 56)
lab = J + np.array([0, 7.0])
ax.text(lab[0], lab[1], r"Phase current $i(t)$: measured" + "\n" + "Module losses: not measured", ha="center", va="bottom",
        color=C["dark"], linespacing=1.2, zorder=12, bbox=dict(fc="white", ec="none", pad=0.8))
ax.plot([lab[0], J[0]], [lab[1] - 0.5, J[1]], color=C["dark"], lw=0.9, zorder=11)
ax.add_patch(Circle(J, 0.9, fc=C["dark"], ec="none", zorder=11))
for mx in (12.0, 88.0):
    ax.add_patch(FancyArrowPatch(J, P(mx + 20, 30.5, PZ + 19), arrowstyle="-|>", mutation_scale=8, lw=0.9,
                                 color=C["dark"], shrinkA=0, shrinkB=1, zorder=11))
# call-outs
s0 = P(88 + 7 + 30, 2.6, PZ); s1 = P(150, -26, 0)
ax.add_patch(FancyArrowPatch(s1, s0, arrowstyle="-|>", mutation_scale=7, lw=0.6, color="black", zorder=12))
ax.text(s1[0], s1[1] - 1.0, "Thermocouples at\nmatched positions", ha="center", va="top", zorder=12)
h0 = P(26, 1, PZ); h1 = P(4, -34, 0)
ax.add_patch(FancyArrowPatch(h1, h0, arrowstyle="-|>", mutation_scale=7, lw=0.6, color=C["red"], zorder=12))
ax.text(h1[0], h1[1] - 1.0, "Heat of module 1 reaches" + chr(10) + "the neighbor only slowly", ha="center", va="top", color=C["red"], zorder=12)
# learned map between matched sensors (the core idea)
m0, m1 = P(104, 1.5, PZ + 0.5), P(40, 1.5, PZ + 0.5)
ax.add_patch(FancyArrowPatch(m0, m1, arrowstyle="-|>", mutation_scale=9, lw=1.2, color=C["purple"], ls=(0, (3, 2)),
                             connectionstyle="arc3,rad=-0.55", zorder=13))
mm = P(76, -44, 0)
ax.text(mm[0], mm[1], r"Learned map $G$:" + chr(10) + "module 2 → module 1", ha="center", va="top", color=C["purple"], linespacing=1.2, zorder=13)
fig.text(0.004, 0.985, "(a)", fontweight="bold", va="top")

# ============================================================ (b) procedure
bx = fig.add_axes([0.445, 0.0, 0.555, 1.0]); bx.set_xlim(0, 100); bx.set_ylim(0, 66.5); bx.axis("off")
fig.text(0.449, 0.985, "(b)", fontweight="bold", va="top")
INK = "#2B2B2B"


def card(x, y, w, h, title, body, col, fc, num=None):
    bx.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.8", fc=fc, ec=col, lw=1.0))
    bx.add_patch(FancyBboxPatch((x, y + h - 6.0), w, 6.0, boxstyle="round,pad=0,rounding_size=1.8", fc=col, ec=col, lw=1.0))
    bx.add_patch(FancyBboxPatch((x, y + h - 6.0), w, 2.5, boxstyle="square,pad=0", fc=col, ec="none"))
    t = (f"{num}  " if num else "") + title
    bx.text(x + w / 2, y + h - 3.1, t, ha="center", va="center", color="white", fontweight="bold", fontsize=9)
    bx.text(x + w / 2, y + (h - 6.0) / 2, body, ha="center", va="center", color=INK, linespacing=1.25)


def arr(p, q, col="#555555"):
    bx.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, lw=1.0, color=col, shrinkA=0, shrinkB=0))


W, H, Y1 = 30.0, 20.0, 40.0
card(1.0, Y1, W, H, "Layout screen", r"$\mathrm{Fo} = \alpha w/\xi$" + "\u00b2" + r"$\,\leq 0.054$" + "\noffline screen, Sec. IV-D", "#6B6B6B", "#F3F3F3", 1)
card(35.0, Y1, W, H, "Baseline", "learn the map\nmodule 2 \u2192 module 1", C["blue"], "#EAF3FA", 2)
card(69.0, Y1, W, H, "Current", "fit one scalar\non the learned map", C["blue"], "#EAF3FA", 3)
card(69.0, 3.0, W, 31.0, "Guards", "G1: response shape\nG2: fast vs wideband\n\nboth two-sided, nominal\nfamily-wise level 5%", "#C98200", "#FFF4E0", 4)
arr((31.0, Y1 + (H - 6) / 2), (35.0, Y1 + (H - 6) / 2)); arr((65.0, Y1 + (H - 6) / 2), (69.0, Y1 + (H - 6) / 2))
arr((84.0, Y1), (84.0, 34.0))
card(35.0, 19.0, W, 15.0, "Pass: report", "drift, nominal 95%\ninterval with allowance", C["green"], "#E6F2EA")
card(35.0, 3.0, W, 15.0, "Fire: abstain", "inconsistency flagged;\ncause not identified", C["red"], "#FBE7DF")
arr((69.0, 26.5), (65.0, 26.5), C["green"]); arr((69.0, 10.5), (65.0, 10.5), C["red"])
# sketch: the drift is the gain change of the learned map
gx = np.linspace(0, 1, 80); x0, y0, wx, hy = 3.0, 3.5, 15.0, 21.0
base = (1 - np.exp(-4 * gx)) / 1.45                                   # baseline response
loss = 1.45 * base                                                    # loss change: same shape, scaled
path = base + 0.30 * np.clip(1 - np.exp(-3.0 * (gx - 0.3)), 0, None) / 1.45   # path change: fast part kept, slow part grows
bx.plot([x0, x0, x0 + wx + 1], [y0 + hy + 1.5, y0, y0], color=INK, lw=0.6)
bx.plot(x0 + wx * gx, y0 + hy * loss, color=C["red"], lw=1.2)
bx.plot(x0 + wx * gx, y0 + hy * path, color="#C98200", lw=1.2, ls=(0, (3, 1.6)))
bx.plot(x0 + wx * gx, y0 + hy * base, color=C["blue"], lw=1.2)
bx.text(x0 + wx + 1.0, y0 + hy * loss[-1] + 0.4, "loss", ha="left", va="center", color=C["red"])
bx.text(x0 + wx + 1.0, y0 + hy * path[-1] - 0.3, "path", ha="left", va="center", color="#C98200")
bx.text(x0 + wx + 1.0, y0 + hy * base[-1] - 1.2, "baseline", ha="left", va="center", color=C["blue"])
bx.text(x0 - 1.0, y0 + hy + 3.0, "loss change: scaled" + chr(10) + "path change: may reshape", ha="left", va="bottom", color=INK, linespacing=1.2)


# fit panel (a) tightly around everything drawn in it (labels included)
for _ in range(4):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    bb = [a.get_window_extent(r) for a in ax.patches + ax.texts if a.get_visible() and a.get_text() != "(a)"
          ] if False else [a.get_window_extent(r) for a in list(ax.patches) + [t for t in ax.texts if t.get_text() != "(a)"]]
    x0 = min(b.x0 for b in bb); x1 = max(b.x1 for b in bb); y0 = min(b.y0 for b in bb); y1 = max(b.y1 for b in bb)
    (dx0, dy0), (dx1, dy1) = inv.transform([(x0, y0), (x1, y1)])
    ax.set_xlim(dx0 - 1.5, dx1 + 1.5); ax.set_ylim(dy0 - 1.5, dy1 + 1.5)
fig.savefig(os.path.join(OUT, "fig_overview_3d.pdf"))
fig.savefig(os.path.join(OUT, "fig_overview_3d.png"), dpi=220)
print("wrote fig_overview_3d")
