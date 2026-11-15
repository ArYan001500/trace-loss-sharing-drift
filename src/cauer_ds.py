"""Datasheet-based thermal generator. Each module is a Cauer ladder obtained from the
manufacturer's Foster junction-to-case network, with a baseplate capacity, the
case-to-sink interface resistance, a heat-sink node and the path to the coolant:

junction ... (datasheet Cauer ladder) ... case node (C_case) [sensor] --R_cs-- sink (C_h) --R_h-- coolant

Loss enters at the junction node and is in watts. Optional first-order sensor lag per reading."""
import numpy as np
from numpy.polynomial import polynomial as Pn
from scipy.linalg import expm


def foster_to_cauer(R, tau):
    """Foster terms (R_i, tau_i) -> Cauer (C_1, R_1, C_2, R_2, ...), junction first."""
    R = np.asarray(R, float); C = np.asarray(tau, float) / R
    num, den = np.zeros(1), np.ones(1)                 # Z = num/den, ascending coefficients in s
    for r, c in zip(R, C):
        num = Pn.polyadd(Pn.polymul(num, [1, r * c]), Pn.polymul(den, [r]))
        den = Pn.polymul(den, [1, r * c])
    num, den = Pn.polytrim(num, 1e-300), Pn.polytrim(den, 1e-300)
    Cs, Rs = [], []
    a, b = den, num                                     # admittance Y = a/b
    for _ in range(len(R)):
        c = a[-1] / b[-1]; Cs.append(c)                 # Y = s c + rem/b
        rem = Pn.polysub(a, Pn.polymul([0, c], b))[:len(b)]
        r = b[-1] / rem[-1]; Rs.append(r)               # b/rem = r + rem2/rem
        rem2 = Pn.polysub(b, Pn.polymul([r], rem))[:len(rem) - 1]
        a, b = rem, rem2
        if len(b) == 0 or np.allclose(b, 0): break
    return np.array(Cs), np.array(Rs)


def chain(Cjc, Rjc, C_case, R_cs, C_h, R_h):
    """State matrices. Nodes: junction ladder C_1 .. C_(n-1), case node C_case (the last Cauer
    capacity is replaced by the baseplate capacity), sink node C_h. The last junction-to-case
    resistance R_n lies in series with the interface R_cs; the sensor reads the baseplate surface,
    between R_n and R_cs. Returns A, B and the output row c (sensor = c @ x)."""
    Cn = np.r_[Cjc[:-1], C_case, C_h]
    Rn = np.r_[Rjc[:-1], Rjc[-1] + R_cs, R_h]
    n = len(Cn); A = np.zeros((n, n)); B = np.zeros((n, 2))
    for i in range(n):
        if i > 0: A[i, i - 1] += 1 / Rn[i - 1]; A[i, i] -= 1 / Rn[i - 1]
        A[i, i] -= 1 / Rn[i]
        if i < n - 1: A[i, i + 1] += 1 / Rn[i]
    B[0, 0] = 1.0; B[-1, 1] = 1 / Rn[-1]
    w = Rjc[-1] / (Rjc[-1] + R_cs); c = np.zeros(n); c[n - 2] = 1 - w; c[n - 1] = w
    return A / Cn[:, None], B / Cn[:, None], c


def discretise(A, B, dt):
    n = A.shape[0]; M = np.zeros((n + 2, n + 2)); M[:n, :n] = A; M[:n, n:] = B
    E = expm(M * dt); return E[:n, :n], E[:n, n:]


def build(p, asym):
    cp = asym.get("cap", 1.0)
    return chain(np.array(p["Cjc"]) * cp, np.array(p["Rjc"]) * asym.get("rda", 1.0), p["C_case"] * cp,
                 p["R_cs"] * asym.get("rb", 1.0), p["C_h"] * cp, p["R_h"] * asym.get("rh", 1.0))


def simulate(I, Tc, p, asym, dt=1.0, alpha=0.006):
    """Two parallel modules sharing phase current I; module 1 carries asym. Returns the
    baseplate-sensor temperatures (n, 2), branch currents and losses (W)."""
    mats = [build(p, asym), build(p, {})]
    sys_ = [discretise(A, B, dt) for A, B, _ in mats]; cs = [m[2] for m in mats]
    s, a, b = p["s"], p["a"], p["b"]
    x = [np.full(m[0].shape[0], Tc[0]) for m in mats]; n = len(I)
    out = np.empty((n, 2)); loss = np.empty((n, 2)); cur = np.empty((n, 2))
    r0 = np.array([1.0 + asym.get("r0", 0.0), 1.0])
    for k in range(n):
        R = r0 * (1 + alpha * (np.array([x[0][0], x[1][0]]) - 25.0)); g = 1 / R; i = I[k] * g / g.sum()
        P = s * (a * (2 * i / 1e3) ** 2 * R + b * (2 * i / 1e3)); P[0] *= asym.get("ploss", 1.0)
        out[k] = [cs[0] @ x[0], cs[1] @ x[1]]; loss[k] = P; cur[k] = i
        for j, (Ad, Bd) in enumerate(sys_):
            x[j] = Ad @ x[j] + Bd @ np.array([P[j], Tc[k]])
    return out, cur, loss


def lag(y, tau, dt=1.0):
    """First-order sensor lag (exact for piecewise-constant input)."""
    if tau <= 0: return y.copy()
    e = np.exp(-dt / tau); z = np.empty_like(y); z[0] = y[0]
    for k in range(1, len(y)): z[k] = e * z[k - 1] + (1 - e) * y[k - 1]
    return z
