%S06_SWEEP_HYSTERESIS  One-parameter sweep in K at Rc = 3.92, upward and downward
%   continuation of the attractor (Finding 6.12, Figure fig_sweep).
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0); P = P0; P.beta = beta_for_Rc(P, 3.92);
Ks = logspace(-2.3, 1.6, 24)*kc; OUT = struct();
for dirn = {'up', 'down'}
    d = dirn{1}; KK = Ks; if strcmp(d, 'down'), KK = fliplr(Ks); end
    res = zeros(numel(KK), 4); y = [];
    for i = 1:numel(KK)
        P.K = KK(i); Xe = upper_ee(P);
        if isempty(y), y = log(Xe.*(1 + 0.05*[1;-1;0;1;1;0])); end
        tt = linspace(400, 700, 6000);
        Yl = ode_long(P, y, tt, 1e-8, 1e-9);
        I = exp(Yl(:,5)); y = Yl(end,:).';
        res(i,:) = [KK(i)/kc, Xe(5), min(I), max(I)];
        fprintf('%s K/Kc=%.4g I*=%.1f min %.3g max %.3g\n', d, res(i,:));
    end
    OUT.(d) = res;
end
save(fullfile(results_dir(), 'sweep.mat'), '-struct', 'OUT');
H = load(fullfile(results_dir(), 'hopf.mat')); kh = H.T(H.T(:,1) == 3.92, [2 5]);
Kg = logspace(-2.5, 1.8, 300); Ie = zeros(size(Kg)); st = false(size(Kg));
for i = 1:numel(Kg), Q = P; Q.K = Kg(i)*kc; [X, ev] = upper_ee(Q); Ie(i) = X(5); st(i) = max(real(ev)) < 0; end
fig = figure('Position', [100 100 560 400]); hold on;
a = Ie; a(~st) = NaN; b = Ie; b(st) = NaN;
loglog(Kg, a, 'k-', 'LineWidth', 1.2); loglog(Kg, b, 'k:');
loglog(OUT.up(:,1), OUT.up(:,4), '^', 'Color', [0.165 0.435 0.592]); loglog(OUT.up(:,1), max(OUT.up(:,3), 1e-12), 'v', 'Color', [0.165 0.435 0.592]);
loglog(OUT.down(:,1), OUT.down(:,4), 'o', 'Color', [0.769 0.271 0.212]); loglog(OUT.down(:,1), max(OUT.down(:,3), 1e-12), 'o', 'Color', [0.769 0.271 0.212]);
for k = kh, plot([k k], [1e-11 5e4], '--', 'Color', [0.6 0.6 0.6]); end
set(gca, 'XScale', 'log', 'YScale', 'log'); ylim([1e-11 5e4]); xlabel('K/K_c'); ylabel('I (per 10^6)'); title('R_c = 3.92');
legend({'I^* stable', 'I^* unstable', 'max I, K increasing', 'min I', 'max I, K decreasing', 'min I'}, 'Location', 'west');
save_fig(fig, 'fig_sweep');
