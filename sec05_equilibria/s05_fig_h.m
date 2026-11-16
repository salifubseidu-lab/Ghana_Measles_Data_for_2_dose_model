%S05_FIG_H  Figure fig_h: the concave function h(lambda) in three regimes.
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0); xs = linspace(0, 0.12, 800);
cases = [0 1.5; 5 0.95; 0.5 0.95]; labs = {'K=0, R_c=1.5','K=5K_c, R_c=0.95','K=0.5K_c, R_c=0.95'};
C = [0.165 0.435 0.592; 0.769 0.271 0.212; 0.416 0.6 0.306];
fig = figure('Position', [100 100 480 360]); hold on;
for j = 1:3
    P = P0; P.K = cases(j,1)*kc; P.beta = beta_for_Rc(P, cases(j,2));
    plot(xs, h_fun(xs, P), 'Color', C(j,:), 'LineWidth', 1.2);
end
for j = 1:3
    P = P0; P.K = cases(j,1)*kc; P.beta = beta_for_Rc(P, cases(j,2)); l = endemic_lambdas(P);
    plot(l, zeros(size(l)), 'o', 'Color', C(j,:), 'MarkerFaceColor', C(j,:));
end
plot(xs, 0*xs, 'k-'); ylim([-0.03 0.04]); xlabel('\lambda (yr^{-1})'); ylabel('h(\lambda)'); legend(labs, 'Location', 'northeast');
save_fig(fig, 'fig_h');
