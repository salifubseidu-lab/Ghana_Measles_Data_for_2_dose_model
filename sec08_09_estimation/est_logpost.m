function lp = est_logpost(M, TH)
%EST_LOGPOST Log posterior (or log-likelihood if M.obj = 'lik', still bounded).
lp = est_logprior(M, TH); ok = isfinite(lp);
if any(ok)
    ll = est_loglik(M, TH(:, ok));
    if strcmp(M.obj, 'lik'), lp(ok) = ll; else, lp(ok) = lp(ok) + ll; end
end
end
