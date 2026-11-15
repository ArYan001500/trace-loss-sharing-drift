function fem_twin(h, hmax, outfile, rb1, src, gap)
% 3D finite-element thermal model of two FMF750DC-66A modules on one aluminium heat sink
% (PDE Toolbox). Independent check of the finite-volume model in src/twin3d.py.
% Layers: aluminium plate 25 mm, grease 1 mm (k chosen so that R'' = 2.7e-4 m^2K/W), baseplate 5 mm (AlSiC).
% Cooling: convection h on the plate bottom (coolant at 0). Loss: 1 W heat flux on the top face of module src.
% Output: rises at the eight sensor points (heat-sink surface under each baseplate, module centre line,
% groove positions 0.69/0.64/0.33/0.28 of the module length), then the same points on the baseplate bottom
% (columns 9-16), sampled every 1 s for 900 s.
if nargin < 1, h = 8000; end
if nargin < 2, hmax = 6e-3; end
if nargin < 3, outfile = sprintf('fem_step_h%d.csv', h); end
if nargin < 4, rb1 = 1; end   % interface-resistance factor of module 1 (grease under module 1)
if nargin < 5, src = 1; end   % module that carries the 1 W step
if nargin < 6, gap = 0.040; end   % distance between the module footprints (m)
mw = 0.086; W = 2*mw + gap + 0.060; L = 0.220; ml = 0.140; y0 = 0.040;
tp = 0.025; tg = 0.001; tb = 0.005;
x1 = (W - 2*mw - gap)/2; x2 = x1 + mw + gap;
R1 = [3 4 0 W W 0 0 0 L L]';
R2 = [3 4 x1 x1+mw x1+mw x1 y0 y0 y0+ml y0+ml]';
R3 = [3 4 x2 x2+mw x2+mw x2 y0 y0 y0+ml y0+ml]';
g2 = decsg([R1 R2 R3], 'R1+R2+R3', char('R1','R2','R3')');
gm = geometryFromEdges(g2);
gm = extrude(gm, tp);
% top faces of the two module regions at z = tp
f1 = nearestFace(gm, [x1+mw/2, y0+ml/2, tp]);
f2 = nearestFace(gm, [x2+mw/2, y0+ml/2, tp]);
gm = extrude(gm, [f1 f2], [tg tb]);
model = createpde('thermal', 'transient');
model.Geometry = gm;
% cell materials by centroid height
nc = gm.NumCells; kc = zeros(nc,1);
% classify cells by height with a coarse probe mesh
pm = createpde('thermal'); pm.Geometry = gm; msh = generateMesh(pm, 'Hmax', 0.01, 'GeometricOrder', 'linear');
for c = 1:nc
    el = findElements(msh, 'region', 'Cell', c);
    zc = mean(msh.Nodes(3, msh.Elements(:, el(1:min(20,numel(el))))), 'all');
    kc(c) = zc; xcen(c) = mean(msh.Nodes(1, msh.Elements(:, el(1:min(20,numel(el))))), 'all'); %#ok<AGROW>
end
for c = 1:nc
    if kc(c) < tp
        thermalProperties(model, 'Cell', c, 'ThermalConductivity', 167, 'MassDensity', 2700, 'SpecificHeat', 900);
    elseif kc(c) < tp + tg
        kg = tg/2.7e-4; if xcen(c) < W/2, kg = kg / rb1; end
        thermalProperties(model, 'Cell', c, 'ThermalConductivity', kg, 'MassDensity', 2500, 'SpecificHeat', 800);
    else
        thermalProperties(model, 'Cell', c, 'ThermalConductivity', 180, 'MassDensity', 3000, 'SpecificHeat', 740);
    end
end
bottom = [];
for f = 1:gm.NumFaces
    n = findNodes(msh, 'region', 'Face', f);
    if all(abs(msh.Nodes(3, n)) < 1e-9), bottom(end+1) = f; end %#ok<AGROW>
end
thermalBC(model, 'Face', bottom, 'ConvectionCoefficient', h, 'AmbientTemperature', 0);
xs = [x1 x2]; ftop1 = nearestFace(gm, [xs(src)+mw/2, y0+ml/2, tp+tg+tb]);
thermalBC(model, 'Face', ftop1, 'HeatFlux', 1/(mw*ml));
thermalIC(model, 0);
generateMesh(model, 'Hmax', hmax, 'GeometricOrder', 'linear');
fprintf('cells %d, bottom faces %d, nodes %d\n', nc, numel(bottom), size(model.Mesh.Nodes, 2));
tlist = 0:1:900;
res = solve(model, tlist);
fr = [0.69 0.64 0.33 0.28];
P = [];
for zs = [tp - 1e-4, tp + tg + 1e-4]          % heat-sink surface (below grease), baseplate bottom (above grease)
    for xc = [x1+mw/2, x2+mw/2]
        for f = fr
            P(:, end+1) = [xc; y0 + f*ml; zs]; %#ok<AGROW>
        end
    end
end
T = interpolateTemperature(res, P(1,:), P(2,:), P(3,:), 1:numel(tlist));
writematrix(T', outfile);
o = (src-1)*4; c = mean(T(5-o:8-o, :), 1) ./ mean(T(1+o:4+o, :), 1);
fprintf('h %d: own rise at 900 s %.2f mK/W; coupling at 4, 10, 30, 120, 900 s: %.4f %.4f %.3f %.3f %.3f\n', ...
    h, 1e3*mean(T(1+o:4+o, end)), c(5), c(11), c(31), c(121), c(901));
end
