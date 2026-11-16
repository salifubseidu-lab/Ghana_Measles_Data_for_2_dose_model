%S10_EQUILIBRIUM_BURDEN  Table 12 (equilibrium burden by coverage, K = 0 and K = 2 Kc)
%   and Figure fig_n_eq (endemic incidence and deaths against p1, K = 0).
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0);
rows = [0.6 0; 0.7 0.6; 0.8 0.6; 0.9 0.6; 0.9 0.9; 0.95 0.6; 0.95 0.9; 0.98 0; 0.98 0.6; 0.98 0.9];
T = zeros(size(rows,1), 9);
fprintf(' p1    p2     Rc   | K=0: cases deaths | K=2Kc: cases deaths stable CFR-share\n');
for r = 1:size(rows,1)
    out = zeros(1,6);
    for kk = 1:2
        P = setcov(P0, rows(r,1), rows(r,2)); P.K = (kk-1)*2*kc; [X, ev] = upper_ee(P);
        if isempty(X), inc = 0; dth = 0; st = 1;
        else, inc = P.eps*X(4); dth = P.d*X(5); st = max(real(ev)) < 0; end
        out(3*kk-2:3*kk) = [inc dth st];
    end
    T(r,:) = [rows(r,:) Rc_fun(P) out(1:2) out(4:6) out(5)/max(out(4),eps)];
    fprintf('%.2f  %.1f  %6.3f | %8.0f %6.0f | %8.0f %6.0f %d  %.3f\n', T(r,:));
end
dlmwrite(fullfile(results_dir(), 'tab12_equilibrium_burden.csv'), T, 'precision', '%.6g');
fig = figure('Position', [100 100 900 340]); p1 = linspace(0.5, 0.999, 300); C = lines(4); p2s = [0 0.3 0.6 0.9];
for j = 1:4
    inc = zeros(size(p1)); dth = inc;
    for i = 1:numel(p1)
        P = setcov(P0, p1(i), p2s(j)); X = upper_ee(P);
        if ~isempty(X), inc(i) = P.eps*X(4); dth(i) = P.d*X(5); end
    end
    subplot(1,2,1); hold on; plot(p1, inc, 'Color', C(j,:));
    subplot(1,2,2); hold on; plot(p1, dth, 'Color', C(j,:));
end
lab = arrayfun(@(v) sprintf('p_2 = %.1f', v), p2s, 'UniformOutput', false);
subplot(1,2,1); xlabel('p_1'); ylabel('cases yr^{-1} per 10^6'); title('(a) endemic incidence (K=0)'); legend(lab);
subplot(1,2,2); xlabel('p_1'); ylabel('deaths yr^{-1} per 10^6'); title('(b) endemic deaths (K=0)'); legend(lab);
save_fig(fig, 'fig_n_eq');
