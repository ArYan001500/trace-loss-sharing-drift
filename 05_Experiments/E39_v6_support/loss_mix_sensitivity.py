"""Post hoc check of assumption A1 (common loss waveform): the transfer-referenced point estimate when the two
modules have different mixes of a current-squared (conduction-like) and a current-proportional (switching-like)
loss term, and the baseline and current periods use different mission phases. The loss laws do not change
between the periods, so the true drift (ratio of mean losses over samples with current above the median) is
close to one; any difference between estimate and truth comes from the violated proportionality.
Finite-element twin step responses of E37, corrected grease (copy of E36; gap 40 mm, h = 1000 W/m2K, sensor below the grease), noise-free.
python loss_mix_sensitivity.py DATA_ROOT"""
import sys, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "src"))
sys.path.insert(0, os.path.join(here, "..", "E05_mission_derivatives"))
from load_mission import load
import tdrift

df = load(sys.argv[1], 0)
PH = {"U": "Ia_ef", "V": "Ib_ef", "W": "Ic_ef"}
for c in PH.values():
    df[c] = df[c].shift(-11)
df = df.dropna(subset=list(PH.values())).reset_index(drop=True)
FEM = os.path.join(here, "..", "E37_fem_twin_grease")
IMP = {s: np.diff(np.loadtxt(os.path.join(FEM, f"fem_g40_h1000_rb1_src{s}.csv"), delimiter=",")[:, :8], axis=0) for s in (1, 2)}


def loss(I, m):
    a = np.abs(I) / 1e3
    q, l = a ** 2, a
    return (1 - m) * q / q.mean() + m * l / l.mean()      # unit-mean mix of the two terms over the mission


def temps(P1, P2):
    out = np.zeros((len(P1), 8))
    for s, P in ((1, P1), (2, P2)):
        x = 450.0 * np.r_[P, P]
        for k in range(8):
            out[:, k] += np.convolve(x, IMP[s][:, k])[len(P):2 * len(P)]
    return 25.0 + out


def truth(I, P1, P2):
    m = np.abs(I) > np.median(np.abs(I))
    return P1[m].mean() / P2[m].mean()


print("baseline  current   m1    m2    true drift   estimate   error %")
for pb, pn in (("U", "V"), ("V", "W"), ("W", "U")):
    Ib, In = df[PH[pb]].to_numpy(), df[PH[pn]].to_numpy()
    for m1, m2 in ((0.0, 0.0), (0.2, 0.0), (0.4, 0.0), (0.4, 0.2), (0.6, 0.0)):
        P1b, P2b, P1n, P2n = loss(Ib, m1), loss(Ib, m2), loss(In, m1), loss(In, m2)
        Tb, Tn = temps(P1b, P2b), temps(P1n, P2n)
        Tc_b, Tc_n = np.full(len(Ib), 25.0), np.full(len(In), 25.0)
        pairs = lambda T: [(f"T{k+1}", T[:, k], T[:, k + 4]) for k in range(4)]
        r = tdrift.estimate_drift(Tc_b, pairs(Tb), Tc_n, pairs(Tn))
        d = truth(In, P1n, P2n) / truth(Ib, P1b, P2b)
        print(f"   {pb}         {pn}     {m1:.1f}   {m2:.1f}    {d:.4f}      {r['drift']:.4f}    {100*(r['drift']/d-1):+.2f}")

# Load-level change: the current-period phase current is scaled by s while the loss laws stay fixed
# (normalization fixed to the baseline mission). A different load level changes the conduction/switching mix.
def loss_fixed(I, m, ref):
    a, ar = np.abs(I) / 1e3, np.abs(ref) / 1e3
    return (1 - m) * a ** 2 / (ar ** 2).mean() + m * a / ar.mean()


print("\nLoad-level change, baseline U, current V scaled by s")
print("  s     m1    m2    true drift   estimate   error %")
Ib, In0 = df[PH["U"]].to_numpy(), df[PH["V"]].to_numpy()
for s in (1.0, 0.6, 0.3):
    In = s * In0
    for m1, m2 in ((0.0, 0.0), (0.2, 0.0), (0.4, 0.0), (0.6, 0.0)):
        P1b, P2b = loss_fixed(Ib, m1, Ib), loss_fixed(Ib, m2, Ib)
        P1n, P2n = loss_fixed(In, m1, Ib), loss_fixed(In, m2, Ib)
        Tb, Tn = temps(P1b, P2b), temps(P1n, P2n)
        pairs = lambda T: [(f"T{k+1}", T[:, k], T[:, k + 4]) for k in range(4)]
        r = tdrift.estimate_drift(np.full(len(Ib), 25.0), pairs(Tb), np.full(len(In), 25.0), pairs(Tn))
        d = truth(In, P1n, P2n) / truth(Ib, P1b, P2b)
        print(f" {s:.1f}   {m1:.1f}   {m2:.1f}    {d:.4f}      {r['drift']:.4f}    {100*(r['drift']/d-1):+.2f}")
