function [x, ok] = newton_ee(x, P, maxit)
%NEWTON_EE Damped Newton iteration for F(x) = 0 using the analytic Jacobian
%   (independent check of the equilibria; no toolbox needed).
if nargin < 3, maxit = 500; end
ok = false;
for it = 1:maxit
    f = model_rhs(x, P);
    if max(abs(f)) < 1e-11*P.Lambda, ok = true; return; end
    dx = -model_jac(x, P)\f; t = 1;
    while any(x + t*dx <= 0) && t > 1e-12, t = t/2; end
    xn = x + t*dx;
    while norm(model_rhs(xn, P)) > norm(f) && t > 1e-10, t = t/2; xn = x + t*dx; end
    x = xn;
end
ok = max(abs(model_rhs(x, P))) < 1e-11*P.Lambda;   % genuine convergence only
end
