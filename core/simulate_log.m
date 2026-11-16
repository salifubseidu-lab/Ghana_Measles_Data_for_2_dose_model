function [t, X] = simulate_log(P, X0, T, n, solver)
%SIMULATE_LOG Integrate model (2.1) in log coordinates on [0,T]; returns X at n points.
if nargin < 4 || isempty(n), n = 4000; end
if nargin < 5, solver = 'ode15s'; end
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9, 'Jacobian', @(t,y) jac_log(t,y,P));
tt = linspace(0, T, n);
f = @(t,y) ode_log(t, y, P);
y0 = log(max(X0(:), 1e-300));
if strcmp(solver, 'ode15s')
    try, [t, Y] = ode15s(f, tt, y0, opts); catch, t = []; end
    if numel(t) < numel(tt) || any(~isfinite(Y(:)))      % failed or incomplete (Octave IDA): fall back
        [t, Y] = ode45(f, tt, y0, odeset('RelTol',1e-8,'AbsTol',1e-9));
    end
else
    [t, Y] = ode45(f, tt, y0, odeset('RelTol',1e-8,'AbsTol',1e-9));
end
X = exp(max(Y, -700)).';      % floor at 1e-304 (exact zero only through underflow)
end
