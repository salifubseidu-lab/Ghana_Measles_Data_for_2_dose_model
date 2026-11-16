function r = RSN_fun(P, tol)
%RSN_FUN Saddle-node threshold R_SN by bisection on the existence of an
%   endemic equilibrium for Rc < 1 (Theorem 5.3).  Returns NaN if K <= Kc.
if nargin < 2, tol = 1e-12; end
if P.K <= Kc_fun(P), r = NaN; return; end
lo = 0.01; hi = 0.999999;
Q = P; Q.beta = beta_for_Rc(P, hi);
if isempty(endemic_lambdas(Q)), r = NaN; return; end
while hi - lo > tol
    mid = (lo + hi)/2; Q.beta = beta_for_Rc(P, mid);
    if ~isempty(endemic_lambdas(Q)), hi = mid; else, lo = mid; end
end
r = hi;
end
