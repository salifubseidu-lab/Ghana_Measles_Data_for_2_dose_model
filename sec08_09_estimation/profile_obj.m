function lp = profile_obj(M, j, val, v)
%PROFILE_OBJ Log-likelihood plus all priors except that of parameter j, at the vector v
%   (all parameters except j) with parameter j fixed at val.
th = [v(1:j-1); val; v(j:end)];
lp = est_logprior(M, th);
if ~isfinite(lp), lp = -1e12; return; end
p = est_unpack(M, th);
switch j                                   % remove the profiled parameter's own prior
    case 1, lp = lp + 0.5*((p.R0 - M.prior_R0(1))/M.prior_R0(2))^2;
    case 5, lp = lp + 0.5*((th(5) - M.RHO_PRIOR(1))/M.RHO_PRIOR(2))^2;
    case 7
        mw = M.prior_omega(1); sw = M.prior_omega(2); vv = mw*(1 - mw)/sw^2 - 1;
        lp = lp - ((mw*vv - 1)*log(p.w1) + ((1 - mw)*vv - 1)*log(1 - p.w1));
    case 9, lp = lp + log(th(9)) + 0.5*((log(th(9)) - log(M.KAPPA_REF))/0.5)^2;
end
lp = lp + est_loglik(M, th);
end
