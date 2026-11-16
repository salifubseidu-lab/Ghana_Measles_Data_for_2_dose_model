function [ll, terms, inc] = est_loglik(M, TH)
%EST_LOGLIK Negative-binomial log-likelihood over observed fitting months (Section 9.3),
%   for each column of TH.  terms: per-month contributions (n_months-by-m).
p = est_unpack(M, TH); inc = est_simulate(M, TH);
mm = inc.*(ones(size(inc,1),1)*p.rho) + 1e-9; k = ones(size(inc,1),1)*p.k; y = M.y*ones(1, size(TH,2));
terms = gammaln(y + k) - gammaln(k) - gammaln(y + 1) + k.*log(k./(k + mm)) + y.*log(mm./(k + mm));
ll = sum(terms(M.mask,:), 1);
ll(~isfinite(ll)) = -1e12;
end
