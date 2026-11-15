"""Comparison of drift estimators on identical scenarios (development seed, not a locked test), frozen method v6.
As E34 compare_baselines.py, with the finite-element twin of E37 (corrected grease) and one added draw per scenario:
with probability 0.5 the module labels are swapped (losses, temperatures and calibration data alike), so that an
interface fault falls on the reference module; the true drift is then inverted.
Generators: 'cauer' (datasheet Cauer model, draws as L3) or 'fem' (finite-element twin, heat-sink sensor,
draws as L5F). Truth: (loss ratio now)/(loss ratio baseline), losses over samples with current above the median.
Methods (all baseline-referenced):
  static    ratio of mean rises r1/r2, now over baseline, median over positions
  wavelet   phase-aware wavelet ratio |S(W1 W2*)/S(|W2|^2)|, periods 20-150 s, median, now over baseline
  twotemp   two-temperature ODE fit (src/twotemp.py, I-13), now over baseline
  slope     slope-jump drift (method v3)
  v6        frozen method v6 (transfer drift with two-sided guards and coupling-aware interval)
  calib     calibrated thermal inversion after Sanz-Alcaine et al. (TPEL 2024; released main.m): linear state model
            T[k+1] = A T[k] + B u[k] identified by least squares from a separate calibration run in which each
            module's loss is stepped independently and MEASURED (labels), then u = pinv(B)(T[k+1] - A T[k]) on
            moving averages (50 s). This baseline uses information the other methods do not have.
python compare_baselines_v6.py ROOT GEN SEED N"""
import sys, os, json
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from cauer_ds import simulate, lag, build, discretise
from load_mission import load
import method_v3, method_v6, twotemp
from wcoh import coherence

root, GEN, SEED, N = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
df = load(root, 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}; CIRC = {"U": "A", "V": "A", "W": "C"}
for c in PH.values(): df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
PAR = json.load(open(os.path.join(here, "..", "E16_datasheet_L3", "ds_params.json")))
FEM = os.path.join(here, "..", "E37_fem_twin_grease")
rng = np.random.default_rng(SEED)

def coolant(ph): return (0.5 * (df[f"Tw_{CIRC[ph]}-In"] + df[f"Tw_{CIRC[ph]}-Out"])).to_numpy()

# ---------- generators: thermal response to given module losses P (n, 2) -> true temperatures (n, 8) -----------
def cauer_temps(P, Tc, lev, asyms):
    out = np.zeros((len(P), 8))
    for j in range(2):
        A, B, c = build(PAR[lev], asyms[j]); Ad, Bd = discretise(A, B, 1.0)
        x = np.full(A.shape[0], Tc[0]); y = np.empty(len(P))
        for k in range(len(P)):
            y[k] = c @ x; x = Ad @ x + Bd @ np.array([P[k, j], Tc[k]])
        out[:, 4 * j:4 * j + 4] = y[:, None]          # one baseplate node per module: same signal at four positions
    return out

IMP = {}
if GEN == "fem":
    for gp in (25, 40, 55):
        for h in (1000, 8000):
            for rb in ("1", "1.3"):
                for s in (1, 2):
                    fn = f"fem_g{gp}_h{h}_rb{rb}_src{s}.csv"
                    IMP[(gp, h, rb, s)] = np.diff(np.loadtxt(os.path.join(FEM, fn), delimiter=",")[:, :8], axis=0)

def fem_temps(P, Tc, gp, h, rb):
    out = np.zeros((len(P), 8))
    for s in (1, 2):
        x = np.r_[P[:, s - 1], P[:, s - 1]]
        for k in range(8):
            out[:, k] += np.convolve(x, IMP[(gp, h, rb, s)][:, k])[len(P):2 * len(P)]
    return Tc[:, None] + out

def readings(T, Tc, g, tl):
    R = np.empty_like(T)
    for k in range(4):
        for j in range(2):
            R[:, k + 4 * j] = lag(Tc + g[k, j] * (T[:, k + 4 * j] - Tc), tl[k, j]) + rng.normal(0, 0.3, len(Tc))
    return R

def pairs(R): return [(f"T{k+1}", R[:, k], R[:, k + 4]) for k in range(4)]

# ---------- estimators ----------
def static(Rb, Tcb, Rn, Tcn):
    q = lambda R, Tc: np.array([np.mean(R[60:, k] - Tc[60:]) / np.mean(R[60:, k + 4] - Tc[60:]) for k in range(4)])
    return float(np.median(q(Rn, Tcn) / q(Rb, Tcb)))

def wavelet(Rb, Tcb, Rn, Tcn):
    def q(R, Tc):
        out = []
        for k in range(4):
            R2, K, P = coherence(R[:, k] - Tc, R[:, k + 4] - Tc)
            band = (P >= 20) & (P <= 150); n = K.shape[1]; e = slice(int(0.15 * n), int(0.85 * n))
            out.append(np.median(np.abs(K[band][:, e])))
        return np.array(out)
    return float(np.median(q(Rn, Tcn) / q(Rb, Tcb)))

def two_temp(Rb, Tcb, Rn, Tcn):
    return twotemp.estimate(Tcn, pairs(Rn))["ratio"] / twotemp.estimate(Tcb, pairs(Rb))["ratio"]

def calib_fit(Rc, Tcc, Pc):
    U = (Rc - Tcc[:, None]).T; X = np.vstack([Pc.T, Tcc[None, :] - Tcc[0]])
    W = U[:, 1:] @ np.linalg.pinv(np.vstack([U[:, :-1], X[:, :-1]]))
    return W[:, :8], W[:, 8:]

def calib_ratio(A, B, R, Tc, I):
    U = (R - Tc[:, None]).T
    mm = lambda Z: np.array([np.convolve(z, np.ones(50) / 50, "same") for z in Z])
    Uf = mm(U)
    est = np.linalg.pinv(B) @ (Uf[:, 1:] - A @ Uf[:, :-1])
    m = (I[:-1] > np.median(I))
    return est[0, m].mean() / est[1, m].mean()

def calib_signal(n, mean_loss):
    P = np.empty((n, 2))
    for j in range(2):
        k = 0
        while k < n:
            L = rng.integers(5, 60); P[k:k + L, j] = mean_loss * rng.uniform(0.2, 1.6); k += L
    return P

# ---------- scenarios ----------
names = ["static", "wavelet", "twotemp", "slope", "v6", "calib"]
rows = []
for s in range(N):
    pb, pn = rng.choice(list("UVW"), 2, replace=False)
    g = rng.normal(1, 0.02, (4, 2)); tl = rng.uniform(0, 3, (4, 2))
    lb, ln = rng.uniform(0.85, 1.15), rng.uniform(0.85, 1.15)
    fault = rng.random() < 0.3
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy(); Tcb, Tcn = coolant(pb), coolant(pn)
    if GEN == "cauer":
        lev = ["0.5", "1.0", "2.0"][rng.integers(3)]; cap = rng.uniform(0.95, 1.05)
        ab = dict(ploss=lb, rh=rng.uniform(0.8, 1.25), cap=cap, rb=1.0)
        an = dict(ploss=ln, rh=rng.uniform(0.8, 1.25), cap=cap, rb=rng.uniform(1.15, 1.40) if fault else 1.0)
        Tb, _, Lb = simulate(Ib, Tcb, PAR[lev], ab); Tn, _, Ln = simulate(In, Tcn, PAR[lev], an)
        Tb = np.repeat(Tb, 4, axis=1); Tn = np.repeat(Tn, 4, axis=1)          # columns: module 1 x4, module 2 x4
        truth = (Ln[In > np.median(In), 0].mean() / Ln[In > np.median(In), 1].mean()) / (Lb[Ib > np.median(Ib), 0].mean() / Lb[Ib > np.median(Ib), 1].mean())
        Pc = calib_signal(1800, Lb.mean()); Tcc = np.full(1800, Tcb.mean())
        Tc_true = cauer_temps(Pc, Tcc, lev, [dict(rh=ab["rh"], cap=cap, rb=1.0), {}])
    else:
        gp = (25, 40, 55)[rng.integers(3)]; h = (1000, 8000)[rng.integers(2)]
        pw = lambda I, l: np.column_stack([l * (I / 1e3) ** 2, (I / 1e3) ** 2]) * 450 / np.mean((I / 1e3) ** 2)
        Pb, Pn = pw(Ib, lb), pw(In, ln)
        Tb = fem_temps(Pb, np.full(len(Ib), 25.0), gp, h, "1"); Tn = fem_temps(Pn, np.full(len(In), 25.0), gp, h, "1.3" if fault else "1")
        Tcb, Tcn = np.full(len(Ib), 25.0), np.full(len(In), 25.0)
        truth = ln / lb
        Pc = calib_signal(1800, 450.0); Tcc = np.full(1800, 25.0)
        Tc_true = fem_temps(Pc, Tcc, gp, h, "1")
    side = 2 if rng.random() < 0.5 else 1
    if side == 2:                                   # swap module labels everywhere
        sw = lambda T: np.concatenate([T[:, 4:8], T[:, 0:4]], axis=1)
        Tb, Tn, Tc_true, Pc, truth = sw(Tb), sw(Tn), sw(Tc_true), Pc[:, ::-1].copy(), 1 / truth
    Rb, Rn = readings(Tb, Tcb, g, tl), readings(Tn, Tcn, g, tl)
    Rc = readings(Tc_true, Tcc, g, tl)
    A, B = calib_fit(Rc, Tcc, Pc)
    v6 = method_v6.estimate_drift((Ib, Tcb, pairs(Rb)), (In, Tcn, pairs(Rn)))
    est = dict(static=static(Rb, Tcb, Rn, Tcn), wavelet=wavelet(Rb, Tcb, Rn, Tcn), twotemp=two_temp(Rb, Tcb, Rn, Tcn),
               slope=v6.get("v3", np.nan), v6=v6.get("drift", np.nan),
               calib=calib_ratio(A, B, Rn, Tcn, In) / calib_ratio(A, B, Rb, Tcb, Ib))
    rows.append((fault, truth, v6["abstain"], side, *[est[k] for k in names]))
    print(f"{s:2d} fault {int(fault)} side {side} truth {truth:.3f} | " + " ".join(f"{k} {100*(est[k]/truth-1):+6.2f}" for k in names)
          + (" | v6 abstains" if v6["abstain"] else ""), flush=True)
R = np.array(rows, dtype=float)
clean, fl, ab = R[:, 0] == 0, R[:, 0] == 1, R[:, 2] == 1
print(f"\n{GEN}: clean n {clean.sum()}, interface n {fl.sum()} (errors in %)")
print(f"{'method':8s} {'clean median':>12s} {'within 3 %':>10s} {'clean max':>9s} {'iface median':>12s}")
for i, k in enumerate(names):
    e = 100 * (R[:, 4 + i] / R[:, 1] - 1)
    print(f"{k:8s} {np.median(np.abs(e[clean])):12.2f} {np.mean(np.abs(e[clean]) <= 3):10.2f} {np.max(np.abs(e[clean])):9.2f} {np.median(e[fl]):+12.1f}")
e = 100 * (R[:, 8] / R[:, 1] - 1)
print(f"v6 reported only: clean n {(clean & ~ab).sum()} median {np.median(np.abs(e[clean & ~ab])):.2f}, within 3 % {np.mean(np.abs(e[clean & ~ab]) <= 3):.2f}; "
      f"interface abstained {ab[fl].mean():.2f}, clean abstained {ab[clean].mean():.2f}")
for sd in (1, 2):
    m = fl & (R[:, 3] == sd)
    print(f"interface faults of module {sd}: n {m.sum()}, v6 abstained {int(ab[m].sum())}; median error % " +
          " ".join(f"{k} {np.median(100 * (R[m, 4 + i] / R[m, 1] - 1)):+.1f}" for i, k in enumerate(names)))
