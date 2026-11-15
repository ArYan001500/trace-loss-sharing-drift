function [drift, ci, logd, v] = tdrift_m(Tcb, B1, B2, Tcn, N1, N2)
% Independent MATLAB implementation of the transfer-referenced drift (point estimate and 95 % interval),
% written from the method description in 03_Protocols/METHOD_v4.md section 3, not translated line by line.
% B1, B2, N1, N2: matrices (samples x positions) of module-1 and module-2 readings, baseline and current.
a = exp(-1/10); M = 6; skip = 60; m = 10;
P = size(B1, 2); logd = zeros(P, 1); v = zeros(P, 1);
for k = 1:P
    b1 = B1(:,k) - Tcb; b2 = B2(:,k) - Tcb; n1 = N1(:,k) - Tcn; n2 = N2(:,k) - Tcn;
    ib = (skip+1:numel(b1))'; in = (skip+1:numel(n1))';
    est = @(ib_, in_) 0.5 * (log(gain(b2, b1, n2, n1, ib_, in_)) - log(gain(b1, b2, n1, n2, ib_, in_)));
    logd(k) = est(ib, in);
    Jb = zeros(m,1); Jn = zeros(m,1);
    cb = chunks(ib, m); cn = chunks(in, m);
    for j = 1:m
        Jb(j) = est(setdiff(ib, cb{j}), in);
        Jn(j) = est(ib, setdiff(in, cn{j}));
    end
    v(k) = max((m-1)/m * sum((Jb - mean(Jb)).^2) + (m-1)/m * sum((Jn - mean(Jn)).^2), 1e-8);
end
% random effects (DerSimonian-Laird tau^2) with Hartung-Knapp-Sidik-Jonkman scale, floor 1
w = 1 ./ v; mu0 = sum(w .* logd) / sum(w); Q = sum(w .* (logd - mu0).^2); df = P - 1;
if df > 0, tau2 = max(0, (Q - df) / (sum(w) - sum(w.^2) / sum(w))); else, tau2 = 0; end
ws = 1 ./ (v + tau2); mu = sum(ws .* logd) / sum(ws);
if P > 1
    q = sum(ws .* (logd - mu).^2) / (P - 1); se = sqrt(max(q, 1) / sum(ws)); crit = tinv(0.975, P - 1);
else
    se = sqrt(1 / sum(ws)); crit = 1.96;
end
drift = exp(mu); ci = exp([mu - crit*se, mu + crit*se]);
end

function g = gain(xb, yb, xn, yn, ib, in)
% learn y = G(q) x + c on the baseline, then fit the scalar on the current period
Xb = lag_basis(xb); th = [Xb(ib,:), ones(numel(ib),1)] \ (yb(ib) - yb(1));
Xn = lag_basis(xn); yh = Xn * th(1:end-1); y = yn - yn(1);
c = [yh(in), ones(numel(in),1)] \ y(in); g = c(1);
end

function X = lag_basis(x)
a = exp(-1/10); M = 6;
x = x - x(1); X = zeros(numel(x), M + 1); X(:,1) = x;
y = filter(sqrt(1 - a^2), [1 -a], x); X(:,2) = y;
for k = 3:M+1
    y = filter([-a 1], [1 -a], y); X(:,k) = y;
end
end

function c = chunks(idx, m)
% contiguous blocks, the first mod(n, m) blocks one element longer
n = numel(idx); base = floor(n/m); extra = mod(n, m); c = cell(m,1); s = 1;
for j = 1:m
    len = base + (j <= extra); c{j} = idx(s:s+len-1); s = s + len;
end
end
