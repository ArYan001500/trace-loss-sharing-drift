"""Builds sanz_1s.npz: 1 s means of the four experiments of the public buck-converter dataset
(github.com/jm-sanz/Automatic-Thermal-Modeling, MIT license; files dataset/Temp_Powers_{Sat,Lin,Drv,Ind}.csv,
14 rows = 8 thermocouple rises and 6 powers, 100 ms sampling). Each output column is the mean of 10 samples.
python make_sanz_1s.py DATASET_DIR OUT.npz"""
import sys, os
import numpy as np
src, out = sys.argv[1], sys.argv[2]
res = {}
for k in ("Sat", "Lin", "Drv", "Ind"):
    X = np.loadtxt(os.path.join(src, f"Temp_Powers_{k}.csv"), delimiter=",")
    n = X.shape[1] // 10
    res[k] = X[:, :10 * n].reshape(X.shape[0], n, 10).mean(axis=2)
np.savez(out, **res)
print({k: v.shape for k, v in res.items()})
