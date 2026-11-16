%S10_HONEYMOON  Finding 10.3 and Figure fig_n_honey: post-vaccination honeymoon after
%   coverage is raised at t = 0 from the baseline equilibrium (K = 0).
setup_paths();
P0 = baseline_params(); X0 = upper_ee(P0);
cov = [0.9 0.6; 0.95 0.6; 0.95 0.9; 0.99 0.9]; C = lines(4);
fig = figure('Position', [100 100 600 380]); hold on; lab = {};
for j = 1:4
    Q = setcov(P0, cov(j,1), cov(j,2)); rc = Rc_fun(Q);
    [t, Y] = simulate_log(Q, X0, 60, 12000); I = Y(5,:);
    semilogy(t, I, 'Color', C(j,:)); lab{end+1} = sprintf('p_1=%.2f, p_2=%.1f (R_c=%.2f)', cov(j,1), cov(j,2), rc); %#ok<SAGROW>
    if rc > 1
        [Imin, imin] = min(I); d = diff(I(imin:end)); k = find(d(1:end-1) > 0 & d(2:end) <= 0, 1) + imin;
        Xe = upper_ee(Q);
        fprintf('(%.2f,%.1f) Rc=%.2f: trough I=%.3g at t=%.2f, first peak at t=%.2f, I*=%.0f\n', cov(j,:), rc, Imin, t(imin), t(k), Xe(5));
    else
        fprintf('(%.2f,%.1f) Rc=%.2f: I(60) = %.3g (elimination)\n', cov(j,:), rc, I(end));
    end
end
plot([0 60], [X0(5) X0(5)], 'k:'); set(gca, 'YScale', 'log'); ylim([1e-12 1e4]);
xlabel('years after coverage change'); ylabel('I per 10^6'); legend(lab, 'Location', 'southwest');
save_fig(fig, 'fig_n_honey');
