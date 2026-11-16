function [lsn, rsn] = fold_point(P)
%FOLD_POINT Fold of the equilibrium branch: minimiser of rho(lambda) (Theorem 6.3).
%   Returns NaN if K <= Kc (rho strictly increasing).
if P.K <= Kc_fun(P), lsn = NaN; rsn = NaN; return; end
lmax = lam_for_rc(P, 1.0);
[lsn, rsn] = fminbnd(@(l) rho_branch(l, P), 1e-12, lmax, optimset('TolX', 1e-14));
end
