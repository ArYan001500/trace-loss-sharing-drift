"""Estimator, method v6 (see 03_Protocols/METHOD_v6.md). Changes from method v5:
V1. Two-sided guards. A path change of either module must be caught, and an interface fault of the reference
    module (module 2) moves both statistics in the direction opposite to a fault of module 1. Both guards form
    one family of two-sided tests at 5 %, each at 2.5 % (Bonferroni): abstain if |z_s| > z_g or |z_x| > z_g,
    z_g = 2.2414 (the 0.9875 quantile of the standard normal distribution).
V2. One sensor set. A position that is unstable over time (|z_stab| > 2 in either period) is excluded from the
    guards and from the point estimate alike.
V3. Guards on the coolant-referenced rise T - Tc, as the point estimate, instead of the raw temperature.
V4. Invalid or non-finite statistics lead to abstention with a reason code; they no longer pass a guard.
V5. Coupling allowance by mapping the interval end points. Settled coupling c between the modules compresses
    the wideband log drift as l_obs = f_c(l) = ln((e^l + c) / (1 + c e^l)). For c in [0, c_max] the true log drift
    of an observed end point y lies between y and f_{c_max}^{-1}(y), with
    f_c^{-1}(y) = ln((e^y - c) / (1 - c e^y)), defined for |y| < ln(1 / c). The reported interval is
    [min(A, f^{-1}(A)), max(B, f^{-1}(B))] for the uncorrected interval [A, B] in logs; c_max = 0.25 as in v5.
    An end point outside the domain of f^{-1} is reported as unbounded (0 or inf).
Point estimate unchanged: transfer-referenced drift (tdrift.py). Slope-jump drift and shape index as method v3.
Inputs per period: phase current I, coolant temperature Tc, and the pairs (name, T1, T2) of module-1 and module-2
readings at identical positions (1 s sampling). Output keys: drift, ci95, ci95_no_coupling, zx, shape_z, abstain,
reason, code (list of reason codes), positions (names used)."""
import numpy as np
from scipy.stats import t as student_t, norm
import method_v3, tdrift

LOCK = dict(z_guard=float(norm.ppf(1 - 0.05 / 4)), c_max=0.25, stab_z=2.0, seed=1)   # load steps: method_v3.LOCK


def _se(ci, k):
    lo, hi = ci
    return np.log(hi / lo) / (2 * (student_t.ppf(0.975, k - 1) if k > 1 else 1.96))


def _uncouple(y, c):
    """Inverse of f_c; returns +-inf outside its domain."""
    e = np.exp(y)
    if e <= c: return -np.inf
    if c * e >= 1: return np.inf
    return float(np.log((e - c) / (1 - c * e)))


def coupling_interval(ci, c_max=None):
    c = LOCK["c_max"] if c_max is None else c_max
    A, B = np.log(ci[0]), np.log(ci[1])
    a, b = min(A, _uncouple(A, c)), max(B, _uncouple(B, c))
    return float(np.exp(a)), float(np.exp(b))


def estimate_drift(base, now):
    """base, now: (I, Tc, pairs) of the baseline and the current period, same sensors."""
    (Ib, Tcb, pb), (In, Tcn, pn) = base, now
    out = dict(flags=[], code=[])

    def stop(code, reason):
        out["code"].append(code); out.update(abstain=True, reason=reason); return out

    rise = lambda Tc, pairs: [(n, a - Tc, b - Tc) for n, a, b in pairs]
    rng = np.random.default_rng(LOCK["seed"])
    zero_b, zero_n = np.zeros(len(Ib)), np.zeros(len(In))
    sb, eb = method_v3._positions(Ib, zero_b, rise(Tcb, pb), rng)
    sn, en = method_v3._positions(In, zero_n, rise(Tcn, pn), rng)
    out["n_steps"] = (len(sb), len(sn))
    if eb is None or en is None:
        return stop("insufficient_excitation", "insufficient excitation (fewer than 20 load steps)")
    used, d, v, f = [], [], [], []
    for name in eb:
        if name not in en: continue
        zb, zn = eb[name]["z_stab"], en[name]["z_stab"]
        if not (np.isfinite(zb) and np.isfinite(zn)) or abs(zb) > LOCK["stab_z"] or abs(zn) > LOCK["stab_z"]:
            out["flags"].append(f"{name}: unstable over time (possible loose sensor), excluded"); continue
        used.append(name)
        d.append(en[name]["log_k"] - eb[name]["log_k"]); v.append(en[name]["v"] + eb[name]["v"])
        f.append((en[name]["shape"] - eb[name]["shape"], en[name]["v_shape"] + eb[name]["v_shape"]))
    out["positions"] = used
    if not used:
        return stop("no_stable_position", "no stable position")
    d, v, fs = np.array(d), np.array(v), np.array(f)
    if not (np.all(np.isfinite(d)) and np.all(np.isfinite(v)) and np.all(np.isfinite(fs)) and np.all(fs[:, 1] > 0)):
        return stop("invalid_statistic", "invalid slope-jump or shape statistic")
    mu, se, tau2, crit = method_v3.random_effects(d, v)
    w = 1 / fs[:, 1]; zs = float((w @ fs[:, 0]) / w.sum() * np.sqrt(w.sum()))
    v3_ci = (np.exp(mu - crit * se), np.exp(mu + crit * se))

    keep = set(used)
    td = tdrift.estimate_drift(Tcb, [p for p in pb if p[0] in keep], Tcn, [p for p in pn if p[0] in keep])
    if not (np.isfinite(td["drift"]) and td["drift"] > 0 and np.all(np.isfinite(td["ci95"])) and td["ci95"][0] > 0):
        return stop("invalid_statistic", "invalid wideband estimate")
    se_td = _se(td["ci95"], len(td["positions"])); se_v3 = _se(v3_ci, len(used))
    zx = float((np.log(td["drift"]) - mu) / np.sqrt(se_td ** 2 + se_v3 ** 2))
    if not (np.isfinite(zs) and np.isfinite(zx)):
        return stop("invalid_statistic", "non-finite guard statistic")

    reasons = []
    if abs(zs) > LOCK["z_guard"]:
        out["code"].append("shape_change"); reasons.append(f"response-shape change near a sensor (shape z {zs:+.2f})")
    if abs(zx) > LOCK["z_guard"]:
        out["code"].append("band_disagreement"); reasons.append(f"fast and wideband drift disagree (z_x {zx:+.2f})")
    out.update(drift=float(td["drift"]), ci95=coupling_interval(td["ci95"]), ci95_no_coupling=tuple(map(float, td["ci95"])),
               v3=float(np.exp(mu)), zx=zx, shape_z=zs, abstain=bool(reasons), reason="; ".join(reasons))
    return out
