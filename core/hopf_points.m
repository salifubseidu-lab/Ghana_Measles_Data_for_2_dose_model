function H = hopf_points(P, rc_hi, n)
%HOPF_POINTS Hopf points on the equilibrium branch for fixed K = P.K, scanning the
%   branch from the fold (or lambda -> 0) up to Rc = rc_hi.
%   Each row of H: [K/Kc, Rc, lambda, l1, omega, I*].
if nargin < 2, rc_hi = 12; end
if nargin < 3, n = 500; end
kc = Kc_fun(P); [lsn, ~] = fold_point(P);
if isnan(lsn), lo = 1e-9; else, lo = lsn*(1 + 1e-7); end
hi = lam_for_rc(P, rc_hi);
lg = logspace(log10(lo), log10(hi), n); sa = arrayfun(@(l) cx_maxre(l, P), lg);
H = zeros(0, 6);
for i = 1:n-1
    if isfinite(sa(i)) && isfinite(sa(i+1)) && sa(i)*sa(i+1) < 0
        lh = fzero(@(l) cx_maxre(l, P), [lg(i) lg(i+1)], optimset('TolX', 1e-16));
        Q = P_on_branch(lh, P); X = ee_from_lambda(lh, Q); [L1, om] = lyap1(X, Q);
        H(end+1,:) = [P.K/kc, rho_branch(lh, P), lh, L1, om, X(5)]; %#ok<AGROW>
    end
end
end
