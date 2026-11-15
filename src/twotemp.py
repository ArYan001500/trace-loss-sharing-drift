"""Two-temperature (input-free) estimator of the loss-sharing ratio of two parallel modules.
Each sensor pair (same position on module 1 and module 2), with r_i = T_i - T_coolant:
    (s + g1)(1 + s tau1) r1 = k (s + g2)(1 + s tau2) r2
fitted through state-variable filters s^j / (s + lam)^3 with bounded parameters.
The loss waveform, the current and the electrical clock are not used."""
import numpy as np
from scipy.signal import bilinear, lfilter
from scipy.optimize import least_squares
from scipy.stats import t as student_t

LOCK = dict(lam=0.08, skip=60, n_jack=10, starts=(0.5, 1.5, 3.0), seed=1)

def _filters(lam):
    den = np.poly([-lam] * 3)
    return [bilinear(np.r_[np.zeros(3 - j), [1.0] + [0.0] * j], den, fs=1.0) for j in range(3)]

def _F(y, j, flt):
    b, a = flt[j]; return lfilter(b, a, y - y[0])

def _fit(A, B, idx):
    def res(q):
        k, t1, t2, g1, g2 = q
        return ((t1 * A[2] + (1 + g1 * t1) * A[1] + g1 * A[0]
                 - k * (t2 * B[2] + (1 + g2 * t2) * B[1] + g2 * B[0])) / (1 + g1 * t1))[idx]
    best = None
    for t0 in LOCK["starts"]:
        r = least_squares(res, [1.0, t0, t0, 0.03, 0.03], bounds=([0.3, 0, 0, 1e-4, 1e-4], [3, 10, 10, 1, 1]))
        if best is None or r.cost < best.cost: best = r
    rss = 2 * best.cost / len(idx)
    return best.x, rss

def position_estimate(r1, r2):
    """Log ratio, its delete-one-block jackknife variance, and the fitted parameters."""
    flt = _filters(LOCK["lam"]); s = LOCK["skip"]
    A = [_F(r1, j, flt)[s:] for j in range(3)]; B = [_F(r2, j, flt)[s:] for j in range(3)]
    n = len(A[0]); full = np.arange(n)
    th, rss = _fit(A, B, full)
    blocks = np.array_split(full, LOCK["n_jack"]); m = len(blocks)
    lj = np.array([np.log(_fit(A, B, np.setdiff1d(full, b))[0][0]) for b in blocks])
    v = (m - 1) / m * np.sum((lj - lj.mean()) ** 2)
    return dict(log_k=np.log(th[0]), v=v, tau1=th[1], tau2=th[2], g1=th[3], g2=th[4], rss=rss)

def random_effects(y, v):
    w = 1 / v; m = (w @ y) / w.sum(); Q = w @ (y - m) ** 2; df = len(y) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if df > 0 else 0.0
    ws = 1 / (v + tau2); mu = (ws @ y) / ws.sum(); k = len(y)
    if k > 1:
        q = (ws @ (y - mu) ** 2) / (k - 1); se = np.sqrt(max(q, 1.0) / ws.sum()); crit = student_t.ppf(0.975, k - 1)
    else:
        se = np.sqrt(1 / ws.sum()); crit = 1.96
    return mu, se, tau2, crit

def estimate(Tc, pairs):
    """pairs: list of (name, T1, T2). Returns the combined ratio, a 95 % interval and per-position fits."""
    pos = {name: position_estimate(T1 - Tc, T2 - Tc) for name, T1, T2 in pairs}
    y = np.array([p["log_k"] for p in pos.values()]); v = np.array([p["v"] for p in pos.values()])
    mu, se, tau2, crit = random_effects(y, v)
    return dict(ratio=np.exp(mu), ci95=(np.exp(mu - crit * se), np.exp(mu + crit * se)), tau2=tau2, positions=pos)
