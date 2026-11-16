function rc = est_rc_year(M, th, years)
%EST_RC_YEAR Control reproduction number Rc(t) at the rates of each calendar year
%   (Theorem 3.3 with iota = 0), for parameter vector th.
p = est_unpack(M, th); D = M.D; rc = zeros(size(years));
for k = 1:numel(years)
    j = find(D.days >= datenum(years(k),1,1), 1);
    mu = D.mu(j); n1 = D.p1(j)*(mu + p.kappa)/(1 - D.p1(j)); n2 = D.p2(j)*(M.ell + mu + p.theta*p.kappa)/(1 - D.p2(j));
    k2 = M.ell + n2 + mu; Dd = mu*(mu + M.ell + n1 + n2) + n1*n2;
    rc(k) = p.R0*mu*(k2 + p.theta*n1)/Dd;
end
end
