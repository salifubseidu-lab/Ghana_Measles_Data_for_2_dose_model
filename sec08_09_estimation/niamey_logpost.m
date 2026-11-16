function lp = niamey_logpost(TH, y, B)
%NIAMEY_LOGPOST Log posterior of the Niamey SEIR fit for each column of TH:
%   NB log-likelihood + uniform prior on rho (logit Jacobian) + R0 ~ N(15, 3^2).
m = size(TH,2); lp = -Inf(1, m);
ok = all(TH >= B(:,1)*ones(1,m) & TH <= B(:,2)*ones(1,m), 1);
if ~any(ok), return; end
T = TH(:,ok); W = numel(y); mm = niamey_sim(T, W) + 1e-9; k = ones(W,1)*exp(T(5,:)); yy = y(:)*ones(1, size(T,2));
ll = gammaln(yy + k) - gammaln(k) - gammaln(yy + 1) + k.*log(k./(k + mm)) + yy.*log(mm./(k + mm));
lp(ok) = sum(ll, 1) + T(3,:) - 2*log1p(exp(T(3,:))) - 0.5*((T(1,:) - 15)/3).^2;
end
