function fem_field_grease(h, hmax, outfile)
% Copy of E28 fem_field.m with the grease area resistance of E37 (R'' = 0.0225 x 0.086 x 0.070 = 1.35e-4 m^2K/W).
% Temperature field on the heat-sink top surface (z = 25 mm - 0.1 mm) after a 1 W loss step in module 1,
% at 4 s, 30 s and 900 s, for figures. Same model as fem_twin.m (gap 40 mm, no interface change).
rb1 = 1; src = 1; gap = 0.040;
% 3D finite-element thermal model of two FMF750DC-66A modules on one aluminium heat sink
% (PDE Toolbox). Independent check of the finite-volume model in src/twin3d.py.
% Layers: aluminium plate 25 mm, grease 1 mm (k chosen so that R'' = 1.35e-4 m^2K/W), baseplate 5 mm (AlSiC).
% Cooling: convection h on the plate bottom (coolant at 0). Loss: 1 W heat flux on the top face of module src.
% Output: rises at the eight sensor points (heat-sink surface under each baseplate, module centre line,
% groove positions 0.69/0.64/0.33/0.28 of the module length), then the same points on the baseplate bottom
% (columns 9-16), sampled every 1 s for 900 s.
if nargin < 1, h = 8000; end
if nargin < 2, hmax = 6e-3; end
if nargin < 3, outfile = sprintf('fem_field_h%d.csv', h); end
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
        kg = tg/(0.0225*0.086*0.070); if xcen(c) < W/2, kg = kg / rb1; end
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
tlist = [0 4 30 900];
res = solve(model, tlist);
[X, Y] = meshgrid(0:0.002:W, 0:0.002:L);
Z = (tp - 1e-4) * ones(size(X));
out = [X(:) Y(:)];
for it = 2:numel(tlist)
    Tq = interpolateTemperature(res, X(:), Y(:), Z(:), it);
    out = [out Tq(:)]; %#ok<AGROW>
end
writematrix(out, outfile);
fprintf('field written: %s (%d points), geometry W %.3f L %.3f x1 %.4f x2 %.4f\n', outfile, size(out,1), W, L, x1, x2);
end
