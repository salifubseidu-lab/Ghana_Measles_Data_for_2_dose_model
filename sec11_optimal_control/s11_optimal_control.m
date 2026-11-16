%S11_OPTIMAL_CONTROL  Table 15, Finding 11.5 and Figures fig_oc_controls/fig_oc_states:
%   optimal strategies none, A (u1), B (u2), C (u3), D (u1,u2), E (u1,u3), F (u2,u3), G (all)
%   for a detected outbreak (I0 = 500), and the early-detection runs (I0 = 10).
%   Set MAXIT (default 3000) to control the optimiser.
setup_paths();
if ~exist('MAXIT', 'var'), MAXIT = 3000; end
P = oc_params(); keys = {'000','100','010','001','110','101','011','111'};
names = {'none','A: u1','B: u2','C: u3','D: u1,u2','E: u1,u3','F: u2,u3','G: all'};
RES = struct();
for scen = {'late', 'early'}
    I0 = 500; if strcmp(scen{1}, 'early'), I0 = 10; end
    X0 = outbreak_ic(P, 1.8, I0, true);
    ks = keys; if strcmp(scen{1}, 'early'), ks = {'000','001','100','111'}; end
    fprintf('\n== %s detection (I0 = %d) ==\n', scen{1}, I0);
    fprintf('%-9s %9s %7s %10s %12s %7s %8s %8s %8s %9s\n', 'strategy','cases','deaths','cost','J','peakI','dose1','dose2','courses','fp-resid');
    for k = 1:numel(ks)
        a = ks{k} == '1'; R = oc_solve(P, X0, a, 'maxit', MAXIT);
        RES.(scen{1}).(['s' ks{k}]) = rmfield(R, 'S');
        fprintf('%-9s %9.0f %7.0f %10.0f %12.0f %7.0f %8.0f %8.0f %8.0f %9.1e\n', names{strcmp(keys, ks{k})}, ...
            R.cases, R.deaths, R.cost, R.J, R.peakI, R.doses1, R.doses2, R.courses, R.fp);
        if any(a)
            on = R.U(:,a) >= 0.999*(ones(size(R.U,1),1)*R.S.umax(a)); off = R.U(:,a) <= 1e-6;
            for j = 1:sum(a)
                lastmax = R.t(find(on(:,j), 1, 'last')); firstoff = R.t(find(~off(:,j), 1, 'last'));
                if isempty(lastmax), lastmax = NaN; end
                fprintf('      control %d: at max until t = %.2f, off by t = %.2f\n', find(cumsum(a) == j, 1), lastmax, firstoff);
            end
        end
    end
end
save(fullfile(results_dir(), 'oc_results.mat'), 'RES');
% figures (late detection)
L = RES.late; C = lines(7); lab = {'u_1','u_2','u_3'};
fig = figure('Position', [100 100 1000 330]);
show = {'s100','s001','s101','s111'}; nm = {'A','C','E','G'};
for j = 1:3
    subplot(1,3,j); hold on;
    for s = 1:numel(show), plot(L.(show{s}).t, L.(show{s}).U(:,j), 'Color', C(s,:)); end
    xlim([0 1.2]); xlabel('time (yr)'); title(lab{j}); legend(nm);
end
save_fig(fig, 'fig_oc_controls');
fig = figure('Position', [100 100 900 340]); show = {'s000','s100','s010','s001','s101','s111'}; nm = {'none','A','B','C','E','G'};
subplot(1,2,1); hold on;
for s = 1:numel(show), semilogy(L.(show{s}).t, L.(show{s}).X(:,5), 'Color', C(s,:)); end
plot([0 5], 1/P.K*[1 1], 'k:'); set(gca, 'YScale', 'log'); xlabel('time (yr)'); ylabel('I'); legend(nm); title('(a) infectious');
subplot(1,2,2); hold on;
for s = 1:numel(show), plot(L.(show{s}).t, cumtrapz(L.(show{s}).t, P.d*L.(show{s}).X(:,5)), 'Color', C(s,:)); end
xlabel('time (yr)'); ylabel('cumulative deaths'); title('(b) deaths');
save_fig(fig, 'fig_oc_states');
