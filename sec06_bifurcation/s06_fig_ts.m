%S06_FIG_TS  Figure fig_ts: (a) damped approach (K = 0) against sustained cycles
%   (K = 1.6 K_H1) at Rc = 3.92; (b) bistability at Rc = 0.9, K = 150 Kc.
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0); H = load(fullfile(results_dir(), 'hopf.mat'));
fig = figure('Position', [100 100 900 340]);
subplot(1,2,1); hold on;
P = P0; P.beta = beta_for_Rc(P, 3.92); x0 = dfe_state(P); X0 = [0.3*x0(1); 0.3*x0(2); x0(3); 50; 50; 0.6e6];
P.K = 0; [t, Y] = simulate_log(P, X0, 150, 8000); semilogy(t, Y(5,:), 'Color', [0.165 0.435 0.592]);
P.K = 1.6*H.KH1; [t, Y] = simulate_log(P, X0, 150, 8000); semilogy(t, Y(5,:), 'Color', [0.769 0.271 0.212]);
set(gca, 'YScale', 'log'); xlabel('time (yr)'); ylabel('I'); legend({'K=0', 'K=1.6K_{H1}'}); title('(a) R_c = 3.92');
subplot(1,2,2); hold on;
P = P0; P.K = 150*kc; P.beta = beta_for_Rc(P, 0.9); x0 = dfe_state(P); Xu = upper_ee(P);
[t, Y] = simulate_log(P, [x0(1:3); 0.2; 0.2; 1.0], 100, 8000); semilogy(t, max(Y(5,:), 1e-30), 'Color', [0.416 0.6 0.306]);
[t, Y] = simulate_log(P, Xu.*[1;1;1;1.2;1.2;1], 100, 8000); semilogy(t, Y(5,:), 'Color', [0.737 0.424 0.145]);
plot([0 100], [Xu(5) Xu(5)], 'k:'); set(gca, 'YScale', 'log'); ylim([1e-12 1e3]); xlabel('time (yr)');
legend({'near DFE, I(0)=0.2', 'near EE, I(0)=1.2I^*'}); title('(b) K=150K_c, R_c = 0.9');
save_fig(fig, 'fig_ts');
