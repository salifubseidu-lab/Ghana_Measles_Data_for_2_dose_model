%S09_SOLVER_IDENTITY_CHECKS  Verification items of Appendix A for Section 9:
%   (1) RK4 (1-day step) against an adaptive solver at the MAP, restarted at every month,
%   pulse and NPI switch; (2) Lemma 8.1 (population closure), (3) Lemma 8.2 (mapping),
%   (4) Proposition 8.3 (pulse map).
setup_paths();
S = load(fullfile(results_dir(), 'est_mcmc.mat')); M = est_setup('covid', true); th = S.R.x; p = est_unpack(M, th); D = M.D;
[inc, st] = est_simulate(M, th);
% (1) adaptive reference solution with drivers held constant within each day
nu1 = D.p1.*(D.mu + p.kappa)./(1 - D.p1); nu2 = D.p2.*(M.ell + D.mu + p.theta*p.kappa)./(1 - D.p2);
k40 = M.K40T + M.DD; beta0 = p.R0*(M.EPS + D.mu(1))*(k40 + D.mu(1))/M.EPS;
x = squeeze(st(1,:,1))'; cut = unique([D.mstart(:); D.sia_idx(:); find(D.t >= 8.25, 1); find(D.t >= 10, 1)]);
sf = min(D.sia_frac*p.h, 0.99); cum = zeros(numel(D.mstart), 1); C = 0;
for c = 1:numel(cut) - 1
    i0 = cut(c); i1 = cut(c+1);
    s = find(D.sia_idx == i0);
    if ~isempty(s), dS = sf(s)*x(1); dV = sf(s)*x(2); x(1) = x(1) - dS; x(2) = x(2) + dS - dV; x(3) = x(3) + dV; end
    rhs = @(t,z) fit_rhs_cont(t, z, D, p, nu1, nu2, beta0, M, k40);
    [~, Z] = ode45(rhs, [D.t(i0) (D.t(i0) + D.t(i1))/2 D.t(i1)], [x; C], odeset('RelTol', 1e-10, 'AbsTol', 1e-6, 'MaxStep', 1/365.25));
    x = Z(end,1:6)'; C = Z(end,7);
    k = find(D.mstart == i1); if ~isempty(k), cum(k) = C; end
end
ref = diff([0; cum(2:end)]); ref = ref(1:numel(inc)); rel = abs(inc - ref)./max(ref, 1);
[~, o] = sort(rel, 'descend');
fprintf('RK4 vs adaptive: median rel. diff %.2e, max %.2e (months %s)\n', median(rel), max(rel), mat2str(o(1:3)'));
% (2) closure
N = squeeze(sum(st(1:end-1,:,1), 2)); P = D.N(D.mstart(1:end-1)); gap = (P - N)./P;
fprintf('(P-N)/P in [%.2e, %.2e]\n', min(gap), max(gap));
% (3) mapping
e1 = max(abs(nu1./(nu1 + p.kappa + D.mu) - D.p1)); e2 = max(abs(nu2./(nu2 + M.ell + p.theta*p.kappa + D.mu) - D.p2));
fprintf('mapping errors: dose 1 %.1e, dose 2 %.1e\n', e1, e2);
% (4) pulse
x = squeeze(st(100,:,1))'; f = 0.5; y = x; y(1) = x(1) - f*x(1); y(2) = x(2) + f*x(1) - f*x(2); y(3) = x(3) + f*x(2);
fprintf('pulse: N error %.1e, dSe %.6g vs formula %.6g\n', abs(sum(y) - sum(x))/sum(x), (x(1) + p.theta*x(2)) - (y(1) + p.theta*y(2)), f*((1 - p.theta)*x(1) + p.theta*x(2)));
