"""Estimator, method v5 (see 03_Protocols/METHOD_v5.md). Identical to method v4 except the guard thresholds:
the two guards form one family of tests at 5 %, so each is tested at 2.5 % (Bonferroni): z < -1.96 and z_x > 1.96.
Interval: thermal coupling between neighbouring modules compresses the wideband drift towards 1 by the factor
(1 - c)/(1 + c) in log, with c the settled coupling (E25 to E28). The 95 % interval is therefore extended on the
side away from 1 by |log d| * 2 c_max / (1 - c_max), c_max = 0.25 (largest settled coupling of the rig twin over
gap >= 25 mm and h >= 1000 W/m2K, E28).
Method v4 header follows.
 Baseline-referenced drift only.
Point estimate and interval: transfer-referenced drift (tdrift.py), which learns the full dynamic map
between the two module rises in the baseline period and fits one scalar on it in the current period.
Guards:
G1. v3 shape index: abstain if the baseline-to-now change of the early/late slope-window index has
    z < -1.96 (v5).
G2. Cross-band consistency: a loss change is frequency-flat, so the wideband drift and the v3 slope-jump
    drift must agree. z_x = (log TD - log v3) / sqrt(se_TD^2 + se_v3^2); abstain if z_x > 1.96 (v5).
    se of each method = log(hi / lo) / (2 t_0.975(k - 1)), k = number of positions used.
G0. Abstain if v3 reports insufficient excitation or no stable position.
Inputs per period: phase current I, coolant temperature Tc, and the pairs (name, T1, T2) of module-1 and
module-2 readings at identical positions (1 s sampling)."""
import numpy as np
from scipy.stats import t as student_t
import method_v3, tdrift

LOCK = dict(zx=1.96, shape_z=-1.96, c_max=0.25)


def _se(ci, k):
    lo, hi = ci
    return np.log(hi / lo) / (2 * (student_t.ppf(0.975, k - 1) if k > 1 else 1.96))


def estimate_drift(base, now):
    """base, now: (I, Tc, pairs)."""
    (Ib, Tcb, pb), (In, Tcn, pn) = base, now
    v3 = method_v3.estimate_drift(base, now)
    out = dict(flags=[], v3=v3.get("drift", np.nan))
    if "drift" not in v3:
        out.update(abstain=True, reason=v3.get("reason", "v3 failed")); return out
    td = tdrift.estimate_drift(Tcb, pb, Tcn, pn)
    se_td = _se(td["ci95"], len(td["positions"])); se_v3 = _se(v3["ci95"], v3["n_used"])
    zx = (np.log(td["drift"]) - np.log(v3["drift"])) / np.sqrt(se_td ** 2 + se_v3 ** 2)
    reasons = []
    if not np.isfinite(v3["shape_z"]) or v3["shape_z"] < LOCK["shape_z"]:
        reasons.append(f"response-shape change near the sensor (shape z {v3['shape_z']:.2f})")
    if zx > LOCK["zx"]:
        reasons.append(f"fast and wideband drift disagree (z_x {zx:.2f})")
    lo, hi = td["ci95"]; ext = abs(np.log(td["drift"])) * 2 * LOCK["c_max"] / (1 - LOCK["c_max"])
    if td["drift"] >= 1: hi = hi * np.exp(ext)
    else: lo = lo * np.exp(-ext)
    out.update(drift=td["drift"], ci95=(lo, hi), ci95_no_coupling=td["ci95"], zx=zx, shape_z=v3["shape_z"], abstain=bool(reasons),
               reason="; ".join(reasons))
    return out
