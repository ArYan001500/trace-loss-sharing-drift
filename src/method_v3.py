"""Estimator, method v3 (see 03_Protocols/METHOD_v3.md). Changes from v2:
C1. The slope jump estimates (loss change / heat capacity). With a junction-to-case ladder from
    a datasheet the capacity tolerance enters one to one (E16), not at 0.6 % as in the fitted
    ladder of E09. The absolute interval therefore carries the full tolerance log(1 + cap_tol).
C2. New baseline-referenced mode (estimate_drift): the ratio now divided by the ratio in a
    baseline period of the same installation. Static nuisances (capacity mismatch, sensor
    coupling, sensor lag mismatch) cancel; the guard is the change of the fast ratio from the
    baseline. The shape index is the log ratio of the module ratios from an early (samples 1-2)
    and a late (samples 4-8) slope window after each step. A pure loss change scales both windows
    alike; a resistance rise in a layer next to the sensor lowers the index. One-sided guard,
    because degradation raises resistance: abstain if the baseline-to-now change of the index is
    below shape_z standard errors (bootstrap over steps, 5 % one-sided). The main ratio keeps the
    window of v2 (samples 1-4).
Inputs per phase: total phase current I, coolant temperature Tc, and for each sensor position
the readings of module 1 and module 2 at identical positions (1 s sampling, aligned clocks).
Output: loss-change ratio of module 1 to module 2 with a 95 % interval, guard flags, and an
abstain decision. No branch current and no cooling state are used."""
import numpy as np
from scipy.stats import t as student_t
from fingerprint_est import load_steps, jumps, fit_two_node

LOCK = dict(min_steps=20, step_frac=0.15, cap_tol=0.05,
            fast_tol=0.06, shape_z=-1.645, main=(1, 5), early=(1, 4), late=(5, 9), stab_z=2.0, n_boot=500, seed=1)

def _ratio(j1, j2):
    return (j1 @ j2) / (j2 @ j2)

def position_estimate(I, Tc, T1, T2, steps, rng):
    j1, j2 = jumps(T1, steps), jumps(T2, steps)
    k = _ratio(j1, j2)
    bs = [ _ratio(j1[s], j2[s]) for s in (rng.integers(0, len(steps), len(steps)) for _ in range(LOCK["n_boot"])) ]
    v = np.var(np.log(bs), ddof=1)
    h = len(steps) // 2                                   # stability: first versus second half of steps
    ka, kb = _ratio(j1[:h], j2[:h]), _ratio(j1[h:], j2[h:])
    z = (np.log(ka) - np.log(kb)) / np.sqrt(2 * v * 2)    # each half has about twice the full-run variance
    p2 = fit_two_node(I, Tc, T2); p1 = fit_two_node(I, Tc, T1, shared=p2[:2])
    k_fast = np.exp(p1[4] - p2[4]); k_slow = np.exp(p1[5] - p2[5])
    k_gain = (p1[2] + p1[3]) / (p2[2] + p2[3])
    return dict(log_k=np.log(k), v=v, z_stab=z, k_fast=k_fast, k_slow=k_slow, k_gain=k_gain)

def random_effects(y, v):
    w = 1 / v; m = (w @ y) / w.sum()
    Q = w @ (y - m) ** 2; df = len(y) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if df > 0 else 0.0
    ws = 1 / (v + tau2); mu = (ws @ y) / ws.sum()
    k = len(y)
    if k > 1:
        q = (ws @ (y - mu) ** 2) / (k - 1)                # HKSJ scale factor
        se = np.sqrt(max(q, 1.0) / ws.sum()); crit = student_t.ppf(0.975, k - 1)
    else:
        se = np.sqrt(1 / ws.sum()); crit = 1.96
    return mu, se, tau2, crit

def _jw(T, steps, a, b):
    d = np.gradient(T)
    return np.array([d[t + a:t + b].mean() - d[t - 5:t - 1].mean() for t in steps])

def drift_position(T1, T2, steps, rng):
    """Main ratio, its bootstrap variance, stability z, and the two-window shape index."""
    n = len(steps); steps = [t for t in steps]
    j1, j2 = _jw(T1, steps, *LOCK["main"]), _jw(T2, steps, *LOCK["main"])
    e1, e2 = _jw(T1, steps, *LOCK["early"]), _jw(T2, steps, *LOCK["early"])
    l1, l2 = _jw(T1, steps, *LOCK["late"]), _jw(T2, steps, *LOCK["late"])
    k = _ratio(j1, j2); sh = np.log(_ratio(e1, e2)) - np.log(_ratio(l1, l2))
    bk, bs = [], []
    for _ in range(LOCK["n_boot"]):
        s = rng.integers(0, n, n)
        bk.append(np.log(_ratio(j1[s], j2[s])))
        re, rl = _ratio(e1[s], e2[s]), _ratio(l1[s], l2[s])
        bs.append(np.log(re) - np.log(rl) if re > 0 and rl > 0 else np.nan)
    v = np.var(bk, ddof=1); vs = np.nanvar(bs, ddof=1)
    h = n // 2; ka, kb = _ratio(j1[:h], j2[:h]), _ratio(j1[h:], j2[h:])
    z = (np.log(ka) - np.log(kb)) / np.sqrt(2 * v * 2)
    return dict(log_k=np.log(k), v=v, z_stab=z, shape=sh, v_shape=vs)

def _positions(I, Tc, pairs, rng):
    steps = load_steps(I, LOCK["step_frac"])
    if len(steps) < LOCK["min_steps"]: return steps, None
    return steps, {name: drift_position(T1, T2, steps, rng) for name, T1, T2 in pairs}

def estimate_drift(base, now):
    """base, now: (I, Tc, pairs) of a baseline period and of the current period, same sensors.
    Output: drift = (loss-change ratio now) / (loss-change ratio in the baseline), 95 % interval."""
    rng = np.random.default_rng(LOCK["seed"]); out = dict(flags=[])
    sb, eb = _positions(*base, rng); sn, en = _positions(*now, rng)
    out["n_steps"] = (len(sb), len(sn))
    if eb is None or en is None:
        out.update(abstain=True, reason="insufficient excitation"); return out
    d, v, f = [], [], []
    for name in eb:
        if abs(eb[name]["z_stab"]) > LOCK["stab_z"] or abs(en[name]["z_stab"]) > LOCK["stab_z"]:
            out["flags"].append(f"{name}: unstable over time (possible loose sensor), excluded"); continue
        d.append(en[name]["log_k"] - eb[name]["log_k"]); v.append(en[name]["v"] + eb[name]["v"])
        f.append((en[name]["shape"] - eb[name]["shape"], en[name]["v_shape"] + eb[name]["v_shape"]))
    if not d:
        out.update(abstain=True, reason="no stable position"); return out
    mu, se, tau2, crit = random_effects(np.array(d), np.array(v))
    fs = np.array(f); w = 1 / fs[:, 1]; ds = (w @ fs[:, 0]) / w.sum(); zs = ds * np.sqrt(w.sum())
    if not np.isfinite(zs) or zs < LOCK["shape_z"]:
        out.update(abstain=True, reason=f"response-shape change near the sensor since baseline (shape z {zs:.2f})")
    else:
        out["abstain"] = False
    half = crit * se
    out.update(drift=np.exp(mu), ci95=(np.exp(mu - half), np.exp(mu + half)), tau2=tau2, shape_change=ds, shape_z=zs, n_used=len(d))
    return out

def estimate(I, Tc, pairs):
    """pairs: list of (name, T1, T2) at identical positions."""
    rng = np.random.default_rng(LOCK["seed"])
    steps = load_steps(I, LOCK["step_frac"])
    out = dict(n_steps=len(steps), positions={}, flags=[])
    if len(steps) < LOCK["min_steps"]:
        out.update(abstain=True, reason="insufficient excitation"); return out
    used = []
    for name, T1, T2 in pairs:
        e = position_estimate(I, Tc, T1, T2, steps, rng); out["positions"][name] = e
        if abs(e["z_stab"]) > LOCK["stab_z"]:
            out["flags"].append(f"{name}: unstable over time (possible loose sensor), excluded"); continue
        used.append(e)
    if not used:
        out.update(abstain=True, reason="no stable position"); return out
    mu, se, tau2, crit = random_effects(np.array([e["log_k"] for e in used]), np.array([e["v"] for e in used]))
    k_fast = np.exp(np.mean([np.log(e["k_fast"]) for e in used]))
    if abs(k_fast - 1) > LOCK["fast_tol"]:
        out.update(abstain=True, reason=f"response-shape change near the sensor (fast ratio {k_fast:.3f})")
    else:
        out["abstain"] = False
    half = crit * se                                      # random part
    sys_ = np.log(1 + LOCK["cap_tol"])                    # heat-capacity tolerance enters one to one (E16)
    out.update(loss_change_ratio=np.exp(mu), ci95=(np.exp(mu - half - sys_), np.exp(mu + half + sys_)), tau2=tau2,
               k_fast=k_fast, k_slow=np.exp(np.mean([np.log(e["k_slow"]) for e in used])),
               k_gain=np.exp(np.mean([np.log(e["k_gain"]) for e in used])), n_used=len(used))
    return out
