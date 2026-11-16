%S09_NIAMEY  Secondary fit (Section 9.11): closed SEIR to the Niamey 2003-04 weekly cases
%   (R0 ~ N(15,3^2), NB likelihood), ensemble MCMC, and the CFR by a Jeffreys binomial posterior.
setup_paths();
if ~exist('NSTEPS_N', 'var'), NSTEPS_N = 12000; end
Dn = read_csv_table(fullfile(data_dir(), 'niamey_weekly.csv')); y = Dn.cases(:); dth = Dn.deaths(:);
B = [3 30; 0.01 0.6; -6 4; log(1) log(500); log(0.5) log(500)];
lp = @(TH) niamey_logpost(TH, y, B);
[x, f] = de_optimize(@(TH) -lp(TH), B, struct('seed', 1, 'maxiter', 400));
[x2, f2] = nm_polish(@(v) -lp(v(:)), x, 5000, 8, 1e-8);     % Nelder-Mead with restarts
if f2 < f, x = x2(:); end
rng(0);
nw = 32; X0 = x*ones(1,nw) + 1e-3*randn(5, nw);
[ch, lpc, acc] = ensemble_mcmc(lp, X0, NSTEPS_N);
burn = floor(NSTEPS_N/3); h = ch(burn+1:10:end,:,:); flat = reshape(permute(h, [1 3 2]), [], 5);
rho = 1./(1 + exp(-flat(:,3))); q = @(v) prctile_(v(:)', [2.5 50 97.5]);
v = q(flat(:,1).*flat(:,2)); fprintf('Re(0) = %.2f [%.2f, %.2f]\n', v(2), v(1), v(3));
v = q(rho); fprintf('reporting = %.3f [%.3f, %.3f]\n', v(2), v(1), v(3));
v = q(rho.*flat(:,2)); fprintf('rho*s0 = %.4f [%.4f, %.4f]\n', v(2), v(1), v(3));
v = q(flat(:,1)); fprintf('R0 (posterior = prior) = %.1f [%.1f, %.1f]\n', v(2), v(1), v(3));
v = q(exp(flat(:,5))); fprintf('NB size = %.2f [%.2f, %.2f];  acceptance %.2f\n', v(2), v(1), v(3), acc);
Dd = sum(dth); Cc = sum(y); cfr = betaincinv([0.025 0.975], Dd + 0.5, Cc - Dd + 0.5);
fprintf('CFR = %.2f%% [%.2f, %.2f] (%d/%d)\n', 100*Dd/Cc, 100*cfr(1), 100*cfr(2), Dd, Cc);
sel = flat(randperm(size(flat,1), 400), :)'; pred = zeros(numel(y), 400);
pred = niamey_sim(sel, numel(y));
pq = prctile_(pred', [2.5 50 97.5])';
fig = figure('Position', [100 100 480 320]); hold on; w = 1:numel(y);
fill([w fliplr(w)], [pq(:,1)' fliplr(pq(:,3)')], [0.8 0.87 0.93], 'EdgeColor', 'none');
plot(w, pq(:,2), 'Color', [0.165 0.435 0.592]); plot(w, y, 'k.', 'MarkerSize', 10);
xlabel('week (1 = week 45 of 2003)'); ylabel('cases per week');
save_fig(fig, 'fig_e_niamey');
save(fullfile(results_dir(), 'est_niamey.mat'), 'flat', 'acc', 'cfr', 'pq');
