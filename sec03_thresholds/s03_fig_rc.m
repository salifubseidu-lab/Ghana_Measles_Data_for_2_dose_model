%S03_FIG_RC  Figure fig_rc: Rc against first-dose uptake p1,
%   (a) for several p2 (theta = 0.15), (b) for several theta (p2 = 0.6).
setup_paths();
P0 = baseline_params(); mu = P0.mu; p1 = linspace(0, 0.999, 400);
C = [0.165 0.435 0.592; 0.769 0.271 0.212; 0.416 0.6 0.306; 0.737 0.424 0.145];
fig = figure('Position', [100 100 900 350]);
subplot(1,2,1); hold on;
p2s = [0 0.3 0.6 0.9];
for j = 1:4
    P = P0; P.nu2 = p2s(j)/(1 - p2s(j))*(P.ell + mu); r = zeros(size(p1));
    for i = 1:numel(p1), P.nu1 = p1(i)/(1 - p1(i))*mu; r(i) = Rc_fun(P); end
    semilogy(p1, r, 'Color', C(j,:), 'LineWidth', 1.2);
end
plot([0 1], [1 1], 'k:'); set(gca, 'YScale', 'log'); xlabel('p_1'); ylabel('R_c');
legend(arrayfun(@(v) sprintf('p_2 = %.1f', v), p2s, 'UniformOutput', false), 'Location', 'southwest'); title('(a) \theta = 0.15');
subplot(1,2,2); hold on;
ths = [0 0.15 0.3 0.5];
for j = 1:4
    P = P0; P.theta = ths(j); r = zeros(size(p1));
    for i = 1:numel(p1), P.nu1 = p1(i)/(1 - p1(i))*mu; r(i) = Rc_fun(P); end
    semilogy(p1, r, 'Color', C(j,:), 'LineWidth', 1.2);
end
plot([0 1], [1 1], 'k:'); set(gca, 'YScale', 'log'); xlabel('p_1');
legend(arrayfun(@(v) sprintf('\\theta = %.2f', v), ths, 'UniformOutput', false), 'Location', 'southwest'); title('(b) p_2 = 0.6');
save_fig(fig, 'fig_rc');
