%S10_BASELINE_TS  Finding 10.1 and Figure fig_n_ts: baseline dynamics from an outbreak
%   initial state (unsaturated effective reproduction number 1.8, I0 = 10) for K = 0
%   and K = 2 Kc: first peak, equilibrium and its eigenvalues, cycle period and the
%   burden averaged over complete limit cycles.
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0);
fig = figure('Position', [100 100 1000 340]); kr = [0 2];
for j = 1:2
    P = P0; P.K = kr(j)*kc; X0 = outbreak_ic(P, 1.8, 10);
    [t, Y] = simulate_log(P, X0, 60, 12001);
    [pk, ip] = max(Y(5,:)); [Xe, ev] = upper_ee(P);
    fprintf('K = %g Kc: first peak I = %.0f at t = %.2f yr; I* = %.0f\n', kr(j), pk, t(ip), Xe(5));
    [~, ii] = max(imag(ev)); fprintf('  leading eigenvalues %.3f +- %.2fi', real(ev(ii)), imag(ev(ii)));
    fprintf(' (period %.2f yr, e-folding %.1f yr)\n', 2*pi/imag(ev(ii)), -1/real(ev(ii)));
    fprintf('  equilibrium burden: %.0f cases/yr, %.0f deaths/yr\n', P.eps*Xe(4), P.d*Xe(5));
    fprintf('  fractions at t=60: S %.3f V1 %.3f V2 %.3f R %.3f\n', Y([1 2 3 6], end)/sum(Y(:,end)));
    if j == 1
        subplot(1,3,1); hold on; idx = [1 2 3 6]; C = lines(4);
        for k = 1:4, plot(t, Y(idx(k),:)./sum(Y,1), 'Color', C(k,:)); end
        legend({'S','V_1','V_2','R'}); xlabel('time (yr)'); title('(a) K=0: fractions');
    end
    subplot(1,3,j+1); hold on;
    semilogy(t, Y(4,:), 'Color', [0.737 0.424 0.145]); semilogy(t, Y(5,:), 'Color', [0.769 0.271 0.212]);
    plot([0 60], [Xe(5) Xe(5)], 'k:'); set(gca, 'YScale', 'log'); ylim([1e-10 5e4]); xlabel('time (yr)');
    legend({'E','I'}); title(sprintf('(%s) K = %g K_c: E, I per 10^6', char('a'+j), kr(j)));
    if j == 2     % burden averaged over complete limit cycles (long run)
        tt = linspace(400, 800, 400001);
        Z = exp(ode_long(P, log(X0), tt, 1e-9, 1e-10)).';
        I = Z(5,:); pkidx = find(I(2:end-1) > I(1:end-2) & I(2:end-1) >= I(3:end) & I(2:end-1) > 1000) + 1;
        a = pkidx(1); b = pkidx(end); per = mean(diff(tt(pkidx)));
        fprintf('  limit cycle: period %.2f yr, max I %.3g, min I %.3g\n', per, max(I), min(I));
        fprintf('  cycle-averaged burden: %.0f cases/yr, %.0f deaths/yr\n', mean(P.eps*Z(4,a:b)), mean(P.d*Z(5,a:b)));
    end
end
save_fig(fig, 'fig_n_ts');
