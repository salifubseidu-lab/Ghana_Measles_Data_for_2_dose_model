function out = basin_outcome(P, X)
%BASIN_OUTCOME 1 if the solution from X ends endemic at t = 150, 0 if I falls below 1e-12.
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'Events', @ext_event, 'Jacobian', @(t,y) jac_log(t,y,P));
[~, y, te] = ode_stiff(@(t,y) ode_log(t,y,P), [0 150], log(X), opts);
if ~isempty(te), out = 0; else, out = double(exp(y(end,5)) > 1e-3); end
end
function [v, term, dir] = ext_event(~, y)
v = y(5) - log(1e-12); term = 1; dir = -1;
end
