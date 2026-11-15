% Post hoc: twin step responses with the grease area resistance derived from the datasheet per half module
% over the full module footprint, R'' = 0.0225 K/W * 0.086 m * 0.070 m = 1.35e-4 m^2K/W (E28 used 2.7e-4).
Rg = 0.0225*0.086*0.070;
for gp = [0.025 0.040 0.055]
  for h = [1000 8000]
    fem_twin_grease(h, 6e-3, sprintf('fem_Rg_g%d_h%d_rb1_src1.csv', round(gp*1e3), h), 1, 1, gp, Rg);
  end
end
for h = [1000 8000]
  fem_twin_grease(h, 6e-3, sprintf('fem_Rg_g40_h%d_rb1.3_src1.csv', h), 1.3, 1, 0.040, Rg);
end
