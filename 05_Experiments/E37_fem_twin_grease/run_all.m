% Finite-element twin with the grease area resistance derived from the datasheet case-to-sink resistance
% per half module (22.5 K/kW) over the full module footprint: R'' = 0.0225 * 0.086 * 0.070 = 1.35e-4 m^2K/W.
% Same geometry, mesh and outputs as E28 fem_twin.m; only the grease resistance differs.
Rg = 0.0225*0.086*0.070;
for gp = [0.025 0.040 0.055]
  for h = [1000 8000]
    for rb = [1 1.3]
      for src = [1 2]
        fn = sprintf('fem_g%d_h%d_rb%s_src%d.csv', round(gp*1e3), h, num2str(rb), src);
        fprintf('%s\n', fn);
        fem_twin_grease(h, 6e-3, fn, rb, src, gp, Rg);
      end
    end
  end
end
