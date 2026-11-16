%S10_VACCINATION_MAPS  Finding 10.2 and Figure fig_n_maps: log10 Rc over (p1,p2) and
%   (p1,omega1); second-dose threshold and first-dose uptake needed for Rc = 1.
setup_paths();
P0 = baseline_params();
nu2min = P0.mu*(P0.theta*R0_fun(P0) - 1);
fprintf('second-dose threshold: nu2 > %.4f, i.e. p2 > %.3f\n', nu2min, nu2min/(nu2min + P0.ell + P0.mu));
for b = [0 0.3 0.6 0.9 0.95]
    f = @(a) Rc_fun(setcov(P0, a, b)) - 1;
    if f(1 - 1e-12) < 0, fprintf('p2 = %.2f: p1* = %.4f\n', b, fzero(f, [0.5 1 - 1e-12]));
    else, fprintf('p2 = %.2f: Rc > 1 for every p1 (Rc(p1=0.98) = %.2f)\n', b, Rc_fun(setcov(P0, 0.98, b))); end
end
p1 = linspace(0.85, 0.9995, 300); p2 = linspace(0, 0.99, 300); w = linspace(0.5, 1, 300);
R1 = zeros(300); R2 = zeros(300);
for i = 1:300
    for j = 1:300
        R1(i,j) = Rc_fun(setcov(P0, p1(j), p2(i)));
        R2(i,j) = Rc_fun(setcov(P0, p1(j), [], w(i)));
    end
end
fig = figure('Position', [100 100 950 380]);
subplot(1,2,1); contourf(p1, p2, log10(R1), linspace(-0.8, 0.8, 17)); hold on;
contour(p1, p2, R1, [1 1], 'k', 'LineWidth', 1.5); contour(p1, p2, R1, [0.5 2 4 8], 'k:'); plot(0.9, 0.6, 'kx');
xlabel('p_1'); ylabel('p_2'); title('(a) R_c over (p_1, p_2)'); caxis([-0.8 0.8]); colorbar;
subplot(1,2,2); contourf(p1, w, log10(R2), linspace(-0.8, 0.8, 17)); hold on;
contour(p1, w, R2, [1 1], 'k', 'LineWidth', 1.5); contour(p1, w, R2, [0.5 2 4 8], 'k:'); plot(0.9, 0.85, 'kx');
xlabel('p_1'); ylabel('\omega_1'); title('(b) R_c over (p_1, \omega_1), p_2 = 0.6'); caxis([-0.8 0.8]); colorbar;
save_fig(fig, 'fig_n_maps');
