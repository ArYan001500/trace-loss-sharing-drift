"""3D finite-volume thermal model of two power modules (FMF750DC-66A footprint, 86 x 140 mm) on one
liquid-cooled aluminium heat sink, as in the traction rig of Zenodo 20395811 (Figs. 4, 5, 7, 8).
Layers, bottom to top: aluminium plate (5 cell layers), grease interface (area resistance), baseplate.
Loss enters the baseplate over the chip area. Sensors: heat-sink top surface under each baseplate, on the
module centre line, at the four groove positions of Fig. 8. Coolant: channels along the plate length at
the two outer edges and between the modules (variant "channels"), or a uniform bottom sink ("bottom").
Backward Euler, sparse LU. Returns step responses of every sensor to a unit (1 W) loss step in module 1."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

GEOM = dict(dx=4e-3, width=272e-3, length=220e-3, mod_w=86e-3, mod_l=140e-3, gap=40e-3, y0=40e-3,
            plate_t=25e-3, n_plate=5, base_t=5e-3, chip_mx=13e-3, chip_my=15e-3,
            k_al=167.0, rc_al=2.43e6, k_b=180.0, rc_b=2.22e6, r_grease=2.7e-4,
            sens_frac=(0.69, 0.64, 0.33, 0.28), chan_w=16e-3, chan_layer=2)


def build(h, cooling="channels", g=None):
    g = dict(GEOM, **(g or {}))
    dx = g["dx"]; nx = int(round(g["width"] / dx)); ny = int(round(g["length"] / dx))
    nzp = g["n_plate"]; dzp = g["plate_t"] / nzp; nz = nzp + 1
    x = (np.arange(nx) + 0.5) * dx; y = (np.arange(ny) + 0.5) * dx
    x1 = (g["width"] - 2 * g["mod_w"] - g["gap"]) / 2
    mods = [(x1, x1 + g["mod_w"]), (x1 + g["mod_w"] + g["gap"], x1 + 2 * g["mod_w"] + g["gap"])]
    y0, y1 = g["y0"], g["y0"] + g["mod_l"]
    inmod = [((x >= a) & (x < b))[:, None] & ((y >= y0) & (y < y1))[None, :] for a, b in mods]
    idx = lambda i, j, k: (k * ny + j) * nx + i
    N = nx * ny * nz
    cap = np.zeros(N); rows, cols, vals = [], [], []
    def link(a, b, G):
        rows.extend([a, b, a, b]); cols.extend([a, b, b, a]); vals.extend([G, G, -G, -G])
    sink = np.zeros(N)
    for k in range(nz):
        base = k == nz - 1
        dz = g["base_t"] if base else dzp
        for j in range(ny):
            for i in range(nx):
                if base and not (inmod[0][i, j] or inmod[1][i, j]):
                    cap[idx(i, j, k)] = 1.0          # dummy cell, isolated
                    continue
                kk, rc = (g["k_b"], g["rc_b"]) if base else (g["k_al"], g["rc_al"])
                cap[idx(i, j, k)] = rc * dx * dx * dz
                if i + 1 < nx and (not base or inmod[0][i + 1, j] or inmod[1][i + 1, j]):
                    link(idx(i, j, k), idx(i + 1, j, k), kk * dz)          # k * (dz*dx) / dx
                if j + 1 < ny and (not base or inmod[0][i, j + 1] or inmod[1][i, j + 1]):
                    link(idx(i, j, k), idx(i, j + 1, k), kk * dz)
                if k + 1 < nzp:
                    link(idx(i, j, k), idx(i, j, k + 1), g["k_al"] * dx * dx / dzp)
                if k == nzp - 1 and (inmod[0][i, j] or inmod[1][i, j]):
                    R = dzp / 2 / (g["k_al"] * dx * dx) + g["r_grease"] / (dx * dx) + g["base_t"] / 2 / (g["k_b"] * dx * dx)
                    link(idx(i, j, k), idx(i, j, k + 1), 1 / R)
    if cooling == "bottom":
        for j in range(ny):
            for i in range(nx):
                sink[idx(i, j, 0)] = h * dx * dx
    else:
        cx = [x1 / 2, x1 + g["mod_w"] + g["gap"] / 2, g["width"] - x1 / 2]
        for c in cx:
            for i in np.where(np.abs(x - c) <= g["chan_w"] / 2)[0]:
                for j in range(ny):
                    sink[idx(i, j, g["chan_layer"])] = h * 2 * dx * dx      # channel top and bottom walls
    K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N)) + sp.diags(sink)
    # loss input: module 1 chip area, top (baseplate) layer
    q = np.zeros(N)
    a, b = mods[0]
    chip = ((x >= a + g["chip_mx"]) & (x < b - g["chip_mx"]))[:, None] & ((y >= y0 + g["chip_my"]) & (y < y1 - g["chip_my"]))[None, :]
    for j in range(ny):
        for i in range(nx):
            if chip[i, j]:
                q[idx(i, j, nz - 1)] = 1.0
    q /= q.sum()
    sens = []
    for m, (a, b) in enumerate(mods):
        i = int(((a + b) / 2) // dx)
        for f in g["sens_frac"]:
            j = int((y0 + f * g["mod_l"]) // dx)
            sens.append(idx(i, j, nzp - 1))
    return dict(K=K, cap=cap, q=q, sens=np.array(sens), N=N)


def step_response(model, T_end=1800.0, dt=0.25):
    """Rise of the 8 sensors (module-1 positions 0..3, module-2 positions 4..7) for a 1 W step in module 1,
    sampled every 1 s."""
    K, C, q = model["K"], model["cap"], model["q"]
    M = (sp.diags(C / dt) + K).tocsc(); lu = splu(M)
    T = np.zeros(model["N"]); out = [np.zeros(len(model["sens"]))]
    per = int(round(1 / dt))
    for n in range(1, int(T_end / dt) + 1):
        T = lu.solve(C / dt * T + q)
        if n % per == 0:
            out.append(T[model["sens"]].copy())
    return np.array(out)
