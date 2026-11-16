%S06_FIG_RHO_BIF  Figures fig_rho (branch function rho(lambda), Theorem 6.3) and
%   fig_bif (equilibrium I* against Rc for K = 0, 150 Kc, 0.5 Kc, coloured by stability).
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0);
C = [0.165 0.435 0.592; 0.769 0.271 0.212; 0.416 0.6 0.306; 0.737 0.424 0.145; 0.482 0.173 0.749];
fig = figure('Position', [100 100 500 370]); hold on; lam = logspace(-6, log10(0.3), 600); krs = [0 0.5 1 5 100];
for j = 1:5, P = P0; P.K = krs(j)*kc; semilogx(lam, rho_branch(lam, P), 'Color', C(j,:), 'LineWidth', 1.2); end
plot([1e-6 0.3], [1 1], ':', 'Color', [0.5 0.5 0.5]); set(gca, 'XScale', 'log'); ylim([0.6 4]);
xlabel('\lambda (yr^{-1})'); ylabel('\rho(\lambda) = R_c on branch');
legend(arrayfun(@(v) sprintf('K = %gK_c', v), krs, 'UniformOutput', false), 'Location', 'northwest');
save_fig(fig, 'fig_rho');
fig = figure('Position', [100 100 1000 340]); krs = [0 150 0.5]; tt = {'(a) K=0', '(b) K=150K_c', '(c) K=0.5K_c'};
for p = 1:3
    subplot(1,3,p); hold on; rcs = linspace(0.5, 4, 700);
    for rc = rcs
        P = P0; P.K = krs(p)*kc; P.beta = beta_for_Rc(P, rc); L = endemic_lambdas(P);
        for k = 1:numel(L)
            X = ee_from_lambda(L(k), P); st = max(real(eig(model_jac(X, P)))) < 0;
            if ~st, col = [0.65 0.65 0.65]; elseif k == numel(L), col = C(1,:); else, col = C(2,:); end
            plot(rc, X(5), '.', 'Color', col, 'MarkerSize', 3);
        end
    end
    plot(rcs(rcs < 1), 0*rcs(rcs < 1), 'k-', 'LineWidth', 1.2); plot(rcs(rcs >= 1), 0*rcs(rcs >= 1), 'k--');
    title(tt{p}); xlabel('R_c'); if p == 1, ylabel('I^* (per 10^6)'); end
    if p == 2, set(gca, 'YScale', 'log'); ylim([1e-2 1e3]); xlim([0.5 2]); end
end
save_fig(fig, 'fig_bif');
