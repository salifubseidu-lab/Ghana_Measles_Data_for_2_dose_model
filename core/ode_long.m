function Y = ode_long(P, y0, tout, rtol, atol)
%ODE_LONG Integrates the log-coordinate model from t = 0 to
%   max(tout) in chunks of at most 50 years and returns the log-state at the times
%   tout (rows of Y).  Uses ode15s with the analytic Jacobian; if the stiff solver
%   fails (e.g. Octave's IDA-based ode15s), the chunk is redone with ode45.
if nargin < 4, rtol = 1e-8; end
if nargin < 5, atol = 1e-9; end
tout = tout(:).'; T = max(tout); edges = unique([0:50:T, T]);
Y = zeros(numel(tout), numel(y0)); y = y0(:);
f = @(t,z) ode_log(t, z, P);
o15 = odeset('RelTol', rtol, 'AbsTol', atol, 'Jacobian', @(t,z) jac_log(t,z,P));
o45 = odeset('RelTol', rtol, 'AbsTol', atol);
for k = 1:numel(edges) - 1
    a = edges(k); b = edges(k+1);
    inside = tout >= a & tout <= b;
    ts = unique([a, tout(inside), b]);
    if numel(ts) == 2, ts = [a, (a + b)/2, b]; end
    try
        [tt, Z] = ode15s(f, ts, y, o15);
        if abs(tt(end) - b) > 1e-9*max(1,b) || any(~isfinite(Z(:))), error('incomplete'); end
    catch
        [tt, Z] = ode45(f, ts, y, o45);
    end
    for j = find(inside), [~, r] = min(abs(tt - tout(j))); Y(j,:) = Z(r,:); end
    y = Z(end,:).';
end
end
