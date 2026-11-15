"""Constant-coefficient fingerprint of two parallel modules from module temperatures,
total phase current and coolant temperature.
k_jump: ratio of slope jumps at load steps; k_gain: ratio of steady gains; k_tau: ratio of
slow time constants (two-node model, loss proxy shared between modules)."""
import numpy as np
from scipy.optimize import least_squares

def _sim(p, I, Tc, x0):
    a, b, k1, k2, lt1, lt2 = p
    P = a * (I / 1e3) ** 2 + b * (I / 1e3); taus = np.exp([lt1, lt2]); ks = np.array([k1, k2])
    x = np.array(x0, float); out = np.empty(len(I))
    for i in range(len(I)):
        out[i] = Tc[i] + x.sum(); x = x + (ks * P[i] - x) / taus
    return out

def fit_two_node(I, Tc, T, shared=None):
    x0 = [T[0] - Tc[0], 0.0]
    if shared is None:
        f = lambda p: _sim(p, I, Tc, x0) - T
        r = least_squares(f, [1, 1, 10, 10, np.log(1), np.log(13)],
                          bounds=([0, 0, 0, 0, np.log(.2), np.log(2)], [50, 50, 200, 200, np.log(5), np.log(500)]))
        return r.x
    a, b = shared
    f = lambda q: _sim([a, b, *q], I, Tc, x0) - T
    r = least_squares(f, [10, 10, np.log(1), np.log(13)],
                      bounds=([0, 0, np.log(.2), np.log(2)], [500, 500, np.log(5), np.log(500)]))
    return np.r_[a, b, r.x]

def load_steps(I, frac=0.15):
    I2 = I ** 2; n = len(I)
    return [t for t in range(6, n - 6) if abs(I2[t + 1] - I2[t - 1]) > frac * I2.max()
            and abs(I2[t + 1] - I2[t - 1]) == max(abs(I2[t + 1 + j] - I2[t - 1 + j]) for j in (-2, -1, 0, 1, 2))]

def jumps(T, steps):
    d = np.gradient(T)
    return np.array([d[t + 1:t + 5].mean() - d[t - 5:t - 1].mean() for t in steps])

def fingerprint(I, Tc, T1, T2):
    st = load_steps(I)
    j1, j2 = jumps(T1, st), jumps(T2, st)
    kj = (j1 @ j2) / (j2 @ j2)
    p2 = fit_two_node(I, Tc, T2); p1 = fit_two_node(I, Tc, T1, shared=p2[:2])
    kg = (p1[2] + p1[3]) / (p2[2] + p2[3])
    kt = np.exp(p1[5]) / np.exp(p2[5])
    return dict(k_jump=kj, k_gain=kg, k_tau=kt, n_steps=len(st), p_ref=p2)
