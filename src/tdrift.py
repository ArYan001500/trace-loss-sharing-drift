"""Transfer-referenced drift estimator.
Baseline period: learn the full dynamic map G from the module-2 rise to the module-1 rise at one
sensor position, r1 = G(q) r2 + c, with an orthonormal Laguerre expansion (pole a) plus a direct term.
G absorbs every static path difference: heat capacity, resistances, sensor coupling and lag.
Current period: fit one scalar d in r1 = d G(q) r2 + c. If the paths are unchanged, d is the change
of the loss ratio since the baseline. Errors-in-variables bias is cancelled to first order by the
geometric mean of the forward fit (r2 -> r1) and the reverse fit (r1 -> r2).
Shape index: d from first differences (fast) minus d from 30 s moving averages (slow), in logs.
A loss change scales all frequencies alike; a resistance rise next to the sensor raises the slow
gain more than the fast gain, so the index falls."""
import numpy as np
from scipy.signal import lfilter
from scipy.stats import t as student_t

LOCK = dict(a=float(np.exp(-1 / 10)), M=6, skip=60, n_jack=10, ma=30, shape_z=-1.645)


def basis(x, a=None, M=None):
    a = LOCK["a"] if a is None else a
    M = LOCK["M"] if M is None else M
    x = x - x[0]
    cols = [x]
    y = lfilter([np.sqrt(1 - a * a)], [1, -a], x)
    cols.append(y)
    for _ in range(M - 1):
        y = lfilter([-a, 1], [1, -a], y)
        cols.append(y)
    return np.column_stack(cols)


def _learn(X, y, idx):
    A = np.column_stack([X[idx], np.ones(len(idx))])
    return np.linalg.lstsq(A, y[idx], rcond=None)[0]


def _gain(yh, y, idx, mode):
    if mode == "lev":
        A = np.column_stack([yh[idx], np.ones(len(idx))])
        return np.linalg.lstsq(A, y[idx], rcond=None)[0][0]
    if mode == "fast":
        u, v = np.diff(yh)[idx[:-1]], np.diff(y)[idx[:-1]]
        return (u @ v) / (u @ u)
    k = np.ones(LOCK["ma"]) / LOCK["ma"]
    u, v = np.convolve(yh, k, "same")[idx], np.convolve(y, k, "same")[idx]
    u, v = u - u.mean(), v - v.mean()
    return (u @ v) / (u @ u)


def _one_dir(xb, yb, xn, yn, ib, i_n):
    th = _learn(basis(xb), yb - yb[0], ib)
    yh = basis(xn) @ th[:-1]
    y = yn - yn[0]
    return {m: _gain(yh, y, i_n, m) for m in ("lev", "fast", "slow")}


def position(b1, b2, n1, n2):
    """b1, b2: baseline rises of module 1 and 2; n1, n2: current rises. Returns the log drift,
    its jackknife variance, the shape index and its variance."""
    s = LOCK["skip"]
    ib, i_n = np.arange(s, len(b1)), np.arange(s, len(n1))

    def est(ib_, in_):
        f = _one_dir(b2, b1, n2, n1, ib_, in_)
        r = _one_dir(b1, b2, n1, n2, ib_, in_)
        ld = 0.5 * (np.log(f["lev"]) - np.log(r["lev"]))
        sh = 0.5 * (np.log(f["fast"]) - np.log(f["slow"]) - np.log(r["fast"]) + np.log(r["slow"]))
        return ld, sh

    ld, sh = est(ib, i_n)
    m = LOCK["n_jack"]
    jb, jn = np.array_split(ib, m), np.array_split(i_n, m)
    Jb = np.array([est(np.setdiff1d(ib, b), i_n) for b in jb])
    Jn = np.array([est(ib, np.setdiff1d(i_n, b)) for b in jn])
    var = lambda J: (m - 1) / m * np.sum((J - J.mean(0)) ** 2, axis=0)
    v = var(Jb) + var(Jn)
    return dict(log_d=ld, v=max(v[0], 1e-8), shape=sh, v_shape=max(v[1], 1e-8))


def random_effects(y, v):
    w = 1 / v
    m = (w @ y) / w.sum()
    Q = w @ (y - m) ** 2
    df = len(y) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if df > 0 else 0.0
    ws = 1 / (v + tau2)
    mu = (ws @ y) / ws.sum()
    k = len(y)
    q = (ws @ (y - mu) ** 2) / (k - 1) if k > 1 else 1.0
    se = np.sqrt(max(q, 1.0) / ws.sum())
    crit = student_t.ppf(0.975, k - 1) if k > 1 else 1.96
    return mu, se, crit


def estimate_drift(Tcb, base_pairs, Tcn, now_pairs):
    """base_pairs / now_pairs: lists of (name, T1, T2) at identical positions."""
    nb = {n: (a - Tcb, b - Tcb) for n, a, b in base_pairs}
    P = {}
    for n, a, b in now_pairs:
        if n in nb:
            P[n] = position(*nb[n], a - Tcn, b - Tcn)
    y = np.array([p["log_d"] for p in P.values()])
    v = np.array([p["v"] for p in P.values()])
    mu, se, crit = random_effects(y, v)
    s = np.array([p["shape"] for p in P.values()])
    vs = np.array([p["v_shape"] for p in P.values()])
    w = 1 / vs
    zs = (w @ s) / w.sum() * np.sqrt(w.sum())
    return dict(drift=np.exp(mu), ci95=(np.exp(mu - crit * se), np.exp(mu + crit * se)), shape_z=zs,
                abstain=bool(zs < LOCK["shape_z"]), positions=P)
