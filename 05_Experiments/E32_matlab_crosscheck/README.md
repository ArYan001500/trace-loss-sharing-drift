# E32 Independent MATLAB implementation of the transfer drift

- `tdrift_m.m` was written from the method description (METHOD_v4.md section 3), not translated line by line. It covers the Laguerre basis, the baseline map, the scalar fit, forward/reverse averaging, the block jackknife and the HKSJ random-effects interval.
- `export_and_compare.py` exports 23 cases and compares the MATLAB results with `src/tdrift.py` on identical CSV inputs:
  - the rig mission halves for U, V and W;
  - 20 Cauer-generator development scenarios (seed 7).
- Result (`crosscheck_output.txt`): the largest |log difference| over the point estimates and both interval ends is 5.2e-13.
