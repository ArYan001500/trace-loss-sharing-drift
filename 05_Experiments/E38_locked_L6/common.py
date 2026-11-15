"""Shared summary of a v6 test run: criteria D1 to D4 of 03_Protocols/METHOD_v6.md and the counts behind D3."""
import numpy as np


def swap(pairs):
    return [(n, b, a) for n, a, b in pairs]


def summary(t, tag):
    clean, fault = ~t.fault, t.fault
    rep = clean & ~t.abstain
    err = t.err_pct.abs()
    D1 = (err[rep] <= 3).mean(); D2 = t.covered[rep].mean()
    D3 = (t.abstain[fault] | (err[fault] <= 3)).mean(); D4 = t.abstain[clean].mean()
    print(f"\n{tag}: clean n {clean.sum()} reported {rep.sum()}; interface n {fault.sum()} "
          f"(module 1: {(fault & (t.side == 1)).sum()}, module 2: {(fault & (t.side == 2)).sum()})")
    for name, v, ok in (("D1 clean |err|<=3 %", D1, D1 >= 0.90), ("D2 coverage", D2, D2 >= 0.85),
                        ("D3 interface abstain or |err|<=3 %", D3, D3 >= 0.80), ("D4 clean abstention", D4, D4 <= 0.10)):
        print(f"{name}: {v:.3f} -> {'PASS' if ok else 'FAIL'}")
    for s in (1, 2):
        m = fault & (t.side == s)
        if m.sum():
            print(f"  interface faults of module {s}: abstained {int(t.abstain[m].sum())}, reported within 3 % "
                  f"{int((~t.abstain[m] & (err[m] <= 3)).sum())}, reported beyond 3 % {int((~t.abstain[m] & (err[m] > 3)).sum())}")
    if rep.sum():
        print(f"clean |err| median {err[rep].median():.2f} %, max {err[rep].max():.2f} %")
    return dict(D1=D1, D2=D2, D3=D3, D4=D4)
