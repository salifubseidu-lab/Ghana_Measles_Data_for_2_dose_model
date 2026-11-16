function [x, f, conv] = nm_polish(fun, x, maxeval, maxrestart, tol)
%NM_POLISH Nelder-Mead (fminsearch) with automatic restarts.
%   Each run gets MAXEVAL evaluations.  The simplex is rebuilt at the current best point
%   and the search restarted until (a) a run ends by its tolerances (exitflag 1) and
%   (b) the restart improves the objective by less than 1e-6, or MAXRESTART is reached.
%   Restarting is the standard remedy for Nelder-Mead stalling or running out of budget
%   in higher dimensions (Lagarias et al. 1998).  conv = true if converged.
if nargin < 3 || isempty(maxeval), maxeval = 5000; end
if nargin < 4 || isempty(maxrestart), maxrestart = 8; end
if nargin < 5 || isempty(tol), tol = 1e-7; end
opts = optimset('MaxFunEvals', maxeval, 'MaxIter', maxeval, 'TolX', tol, 'TolFun', tol, 'Display', 'off');
[x, f, flag] = fminsearch(fun, x, opts); conv = false;
for r = 1:maxrestart
    [x2, f2, flag2] = fminsearch(fun, x, opts);
    gain = f - f2;
    if f2 < f, x = x2; f = f2; end
    if flag2 == 1 && gain < 1e-6, conv = true; break; end
    flag = flag2;
end
if ~conv
    warning('nm_polish:notConverged', ['Nelder-Mead did not converge after %d restarts ', ...
            '(objective %.6f); increase MAXEVAL or MAXRESTART.'], maxrestart, f);
end
end
