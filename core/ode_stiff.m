function [t, y, te, ye, ie] = ode_stiff(f, tspan, y0, opts)
%ODE_STIFF ode15s with an automatic ode45 fallback.  In MATLAB, ode15s is used; the
%   fallback is only for Octave, whose IDA-based ode15s sometimes stops at t = 0 on this
%   model.  Event outputs (te, ye, ie) are requested only when an Events function is set,
%   because MATLAB does not assign them otherwise.
te = []; ye = []; ie = [];
hasEv = ~isempty(odeget(opts, 'Events'));
ok = false;
try
    if hasEv, [t, y, te, ye, ie] = ode15s(f, tspan, y0, opts);
    else,     [t, y] = ode15s(f, tspan, y0, opts); end
    ok = ~isempty(t) && all(isfinite(y(:))) && ...
         (abs(t(end) - tspan(end)) < 1e-9*max(1, abs(tspan(end))) || ~isempty(te));
catch
end
if ~ok
    o = odeset(opts, 'Jacobian', []);
    if hasEv, [t, y, te, ye, ie] = ode45(f, tspan, y0, o);
    else,     [t, y] = ode45(f, tspan, y0, o); end
end
end
