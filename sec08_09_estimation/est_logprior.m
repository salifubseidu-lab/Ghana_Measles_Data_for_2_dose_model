function lp = est_logprior(M, TH)
%EST_LOGPRIOR Log prior (up to a constant) of Table 13; -Inf outside the search bounds.
m = size(TH,2); lp = zeros(1, m); p = est_unpack(M, TH);
out = any(TH < M.bounds(:,1)*ones(1,m) | TH > M.bounds(:,2)*ones(1,m), 1);
lp = lp - 0.5*((p.R0 - M.prior_R0(1))/M.prior_R0(2)).^2;
mw = M.prior_omega(1); sw = M.prior_omega(2); v = mw*(1 - mw)/sw^2 - 1; A = mw*v; B = (1 - mw)*v;
lp = lp + (A - 1)*log(p.w1) + (B - 1)*log(1 - p.w1);
lp = lp - 0.5*((TH(5,:) - M.RHO_PRIOR(1))/M.RHO_PRIOR(2)).^2;
if isempty(M.kappa_fixed)
    kap = TH(9,:); bad = kap <= 0; kap(bad) = 1;
    lp = lp - log(kap) - 0.5*((log(kap) - log(M.KAPPA_REF))/0.5).^2; lp(bad) = -Inf;
end
lp(out) = -Inf;
end
