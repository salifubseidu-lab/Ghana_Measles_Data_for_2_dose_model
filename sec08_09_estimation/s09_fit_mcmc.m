%S09_FIT_MCMC  Section 9.5-9.6: MAP by differential evolution (3 seeds) + Nelder-Mead,
%   AIC comparison with the no-NPI model, and the posterior by ensemble MCMC
%   (44 walkers x NSTEPS iterations; first half discarded).  Diagnostics: acceptance,
%   integrated autocorrelation time, split-Rhat, effective sample size.
setup_paths();
if ~exist('NSTEPS', 'var'), NSTEPS = 10000; end
if ~exist('DEITER', 'var'), DEITER = 600; end
if ~exist('SEEDS', 'var'), SEEDS = [1 2 3]; end
if ~exist('NMEVAL', 'var'), NMEVAL = 5000; end      % per Nelder-Mead run (restarted)
M = est_setup('covid', true);
fprintf('Reference model (NPI term, SIA reach):\n'); R = est_fit(M, SEEDS, DEITER, NMEVAL);
fprintf('MAP loglik %.2f, AIC %.1f\n', R.loglik, R.AIC);
fprintf('No-NPI model:\n'); M0 = est_setup('covid', false); R0 = est_fit(M0, SEEDS(1:min(2,end)), DEITER, NMEVAL);
fprintf('no-NPI loglik %.2f, AIC %.1f, dAIC = %.1f\n', R0.loglik, R0.AIC, R0.AIC - R.AIC);
% MCMC
d = M.npar; nw = 4*d; rng(7);
X0 = R.x*ones(1,nw) + 1e-4*randn(d,nw).*(max(abs(R.x),0.1)*ones(1,nw));
X0 = min(max(X0, (M.bounds(:,1) + 1e-9)*ones(1,nw)), (M.bounds(:,2) - 1e-9)*ones(1,nw));
tic; [chain, lp, acc] = ensemble_mcmc(@(TH) est_logpost(M, TH), X0, NSTEPS); fprintf('MCMC time %.0f s\n', toc);
burn = floor(NSTEPS/2); h = chain(burn+1:end,:,:);
tau = zeros(1,d); for j = 1:d, tau(j) = iat(squeeze(h(:,j,:))); end
n2 = floor(size(h,1)/2); seg = cat(3, h(1:n2,:,:), h(n2+1:2*n2,:,:));    % split chains
W = squeeze(mean(var(seg, 0, 1), 3)); Bv = n2*var(squeeze(mean(seg, 1)), 0, 2)';
rhat = sqrt(((n2 - 1)/n2*W + Bv/n2)./W);
thin = max(1, floor(min(tau)/2)); flat = reshape(permute(h(1:thin:end,:,:), [1 3 2]), [], d);
ess = size(h,1)*nw/max(tau);
fprintf('acceptance %.3f, tau %.1f-%.1f, max split-Rhat %.3f, ESS %.0f, retained %d\n', acc, min(tau), max(tau), max(rhat), ess, size(flat,1));
save(fullfile(results_dir(), 'est_mcmc.mat'), 'R', 'R0', 'flat', 'lp', 'acc', 'tau', 'rhat', 'ess', 'burn', 'nw', 'NSTEPS');
