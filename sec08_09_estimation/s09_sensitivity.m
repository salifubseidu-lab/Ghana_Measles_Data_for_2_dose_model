%S09_SENSITIVITY  Table 20 and Finding 9.6: MAP under alternative data choices and fixed
%   values; Rc(2024) at 2024 rates.
setup_paths();
if ~exist('DEITER', 'var'), DEITER = 600; end
if ~exist('NMEVAL', 'var'), NMEVAL = 5000; end      % per Nelder-Mead run (restarted)
sc = {'Reference (with NPI)', {};  'h_SIA = 1 fixed', {'sia_est', false};  'No NPI term', {'covid', false};
      'Lab+epi-linked only', {'target', 'conf'};  'DHS-anchored coverage', {'cov', 'dhs'};
      'ell = 0', {'ell', 0};  'ell = 0.03', {'ell', 0.03};  'lambda_ref = 1/5.2 fixed', {'kappa_fixed', 1/5.2};
      'flat R0 prior', {'prior_R0', [15 1000]};  'wide omega1 prior', {'prior_omega', [0.75 0.15]};  'MLE (no priors)', {'obj', 'lik'}};
T = zeros(size(sc,1), 9);
for s = 1:size(sc,1)
    args = sc{s,2}; if ~any(strcmp(args(1:2:end), 'covid')), args = [args, {'covid', true}]; end %#ok<AGROW>
    M = est_setup(args{:}); fprintf('%s\n', sc{s,1}); R = est_fit(M, [1 2], DEITER, NMEVAL);
    p = est_unpack(M, R.x);
    T(s,:) = [p.R0 p.w1 p.rho p.kappa p.s0 p.h est_rc_year(M, R.x, 2024) R.loglik R.AIC];
    fprintf('  R0 %.2f  omega1 %.3f  rho %.3f  lambda_ref %.4f  s0 %.4f  h %.3f  Rc(2024) %.3f  loglik %.1f\n', T(s,1:8));
end
save(fullfile(results_dir(), 'est_sensitivity.mat'), 'T', 'sc');
