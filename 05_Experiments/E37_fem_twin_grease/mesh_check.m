% Mesh check of the corrected twin: h 8000, gap 40 mm, loss in module 1, default (6 mm) and finer (4 mm) mesh.
Rg = 0.0225*0.086*0.070;
fem_twin_grease(8000, 6e-3, 'mesh_h6.csv', 1, 1, 0.040, Rg);
fem_twin_grease(8000, 4e-3, 'mesh_h4.csv', 1, 1, 0.040, Rg);
a = readmatrix('mesh_h6.csv'); b = readmatrix('mesh_h4.csv');
fprintf('own rise at 900 s: %.4f and %.4f mK/W (difference %.2f %%)\n', 1e3*mean(a(end,1:4)), 1e3*mean(b(end,1:4)), 100*(mean(b(end,1:4))/mean(a(end,1:4))-1));
fprintf('settled coupling: %.4f and %.4f\n', mean(a(end,5:8))/mean(a(end,1:4)), mean(b(end,5:8))/mean(b(end,1:4)));
fprintf('largest relative difference of the own rise over 1-900 s: %.2f %%\n', 100*max(abs(mean(b(2:end,1:4),2)./mean(a(2:end,1:4),2)-1)));
