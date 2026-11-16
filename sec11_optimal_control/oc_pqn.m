function [z, info] = oc_pqn(fun, z, maxit, m, tol)
%OC_PQN Projected limited-memory quasi-Newton method for min f(z), 0 <= z <= 1
%   (no toolbox needed).  Direction: L-BFGS two-loop recursion on the free variables;
%   step: projected backtracking (Armijo) line search.
if nargin < 3, maxit = 3000; end
if nargin < 4, m = 30; end
if nargin < 5, tol = 1e-10; end
proj = @(v) min(max(v, 0), 1);
[f, g] = fun(z); Sm = zeros(numel(z), 0); Ym = Sm; info.nf = 1;
for it = 1:maxit
    pg = z - proj(z - g);
    if max(abs(pg)) < tol, break; end
    free = ~((z <= 0 & g > 0) | (z >= 1 & g < 0));
    q = g.*free; a = zeros(1, size(Sm,2)); rho = a;
    for i = size(Sm,2):-1:1
        rho(i) = 1/(Ym(:,i)'*Sm(:,i)); a(i) = rho(i)*(Sm(:,i)'*q); q = q - a(i)*Ym(:,i).*free;
    end
    if ~isempty(Sm), gam = (Sm(:,end)'*Ym(:,end))/(Ym(:,end)'*Ym(:,end)); else, gam = 1/max(norm(g.*free, inf), 1e-12); end
    r = gam*q;
    for i = 1:size(Sm,2)
        b = rho(i)*(Ym(:,i)'*r); r = r + (a(i) - b)*Sm(:,i).*free;
    end
    d = -r.*free;
    if g'*d >= 0, d = -g.*free; Sm = Sm(:,[]); Ym = Ym(:,[]); end
    step = 1; accepted = false;
    for ls = 1:40
        zn = proj(z + step*d); [fn, gn] = fun(zn); info.nf = info.nf + 1;
        if fn <= f + 1e-4*(g'*(zn - z)), accepted = true; break; end
        step = step/2;
    end
    if ~accepted, break; end
    s = zn - z; y = gn - g;
    if s'*y > 1e-12*norm(s)*norm(y)
        Sm = [Sm s]; Ym = [Ym y]; %#ok<AGROW>
        if size(Sm,2) > m, Sm(:,1) = []; Ym(:,1) = []; end
    end
    if abs(f - fn) <= 1e-15*max(1, abs(f)) && max(abs(s)) < 1e-12, z = zn; f = fn; g = gn; break; end
    z = zn; f = fn; g = gn;
end
info.it = it; info.f = f; info.pg = max(abs(z - proj(z - g)));
end
