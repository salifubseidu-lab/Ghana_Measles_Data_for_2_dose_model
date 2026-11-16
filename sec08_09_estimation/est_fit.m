function R = est_fit(M, seeds, maxiter, nmeval)
%EST_FIT MAP (or MLE if M.obj = 'lik') by differential evolution from several seeds,
%   each polished by Nelder-Mead with restarts (NM_POLISH).  Returns the best solution, its
%   log-likelihood, AIC and the per-seed objective values.
if nargin < 2, seeds = [1 2 3]; end
if nargin < 3, maxiter = 600; end
if nargin < 4, nmeval = 5000; end          % evaluations per Nelder-Mead run
negobj = @(TH) -est_logpost(M, TH);
R.runs = zeros(numel(seeds), 1); best = Inf;
for s = 1:numel(seeds)
    [x, fx] = de_optimize(negobj, M.bounds, struct('seed', seeds(s), 'maxiter', maxiter));
    pen = @(v) negobj(v(:));
    [x2, f2, conv] = nm_polish(pen, x, nmeval, 8, 1e-7);      % Nelder-Mead with restarts
    if f2 < fx, x = x2(:); fx = f2; end
    R.converged(s) = conv;
    R.runs(s) = fx; fprintf('  seed %d: -objective %.4f (Nelder-Mead converged: %d)\n', seeds(s), fx, conv);
    if fx < best, best = fx; R.x = x; end
end
R.negobj = best;
fprintf('  spread of the objective across seeds: %.2e\n', max(R.runs) - min(R.runs)); R.loglik = est_loglik(M, R.x); R.npar = M.npar; R.AIC = 2*M.npar - 2*R.loglik;
R.nobs = sum(M.mask); R.names = M.names;
end
