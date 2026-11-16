%S09_POSTERIOR_ANALYSIS  Tables 14, 17-19 and Findings 9.1-9.4: posterior summaries,
%   posterior predictive checks (fit and hold-out), CRPS, seasonal-persistence benchmark,
%   annual Rc(t), mean R_eff, infections, and SIA counterfactuals (f_j = 0).
setup_paths();
if ~exist('NDRAW', 'var'), NDRAW = 1500; end
if ~exist('NCF', 'var'), NCF = 300; end
S = load(fullfile(results_dir(), 'est_mcmc.mat')); M = est_setup('covid', true);
rng(11);
flat = S.flat; idx = randperm(size(flat,1), min(NDRAW, size(flat,1))); TH = flat(idx,:)';
q = @(v) prctile_(v, [2.5 50 97.5]);
p = est_unpack(M, TH);
nat = struct('R0',p.R0,'a',p.a,'iota',p.iota,'rho',p.rho,'s0',p.s0,'omega1',p.w1,'theta',p.theta, ...
             'k',p.k,'kappa',p.kappa,'c',p.c,'h',p.h,'peak_month',mod(p.phi,1)*12 + 1);
fn = fieldnames(nat); fprintf('Posterior medians [95%% CrI]:\n');
for i = 1:numel(fn), v = q(nat.(fn{i})); fprintf('  %-10s %.4g [%.4g, %.4g]\n', fn{i}, v(2), v(1), v(3)); end
[inc, st] = est_simulate(M, TH); nm = size(inc,1);
ppc = nb_rand(ones(nm,1)*p.k, inc.*(ones(nm,1)*p.rho));
y = M.D.y; obs = ~isnan(y); fitm = M.D.fit & obs; hold = ~M.D.fit & obs;
lo = prctile_(ppc', 2.5)'; hi = prctile_(ppc', 97.5)'; lo50 = prctile_(ppc', 25)'; hi50 = prctile_(ppc', 75)'; med = prctile_(ppc', 50)';
in95 = y >= lo & y <= hi; in50 = y >= lo50 & y <= hi50;
cr = nan(nm,1); for i = find(obs)', cr(i) = crps_sample(ppc(i, 1:min(400,end)), y(i)); end
naive = [nan(12,1); y(1:end-12)]; hm = hold & ~isnan(naive);
fprintf('\nFit window:  n=%d  cov95 %.1f%%  cov50 %.1f%%  MAE %.1f  CRPS %.1f\n', sum(fitm), 100*mean(in95(fitm)), 100*mean(in50(fitm)), mean(abs(med(fitm) - y(fitm))), mean(cr(fitm)));
fprintf('Hold-out:    n=%d  cov95 %.1f%%  cov50 %.1f%%  MAE %.1f  CRPS %.1f\n', sum(hold), 100*mean(in95(hold)), 100*mean(in50(hold)), mean(abs(med(hold) - y(hold))), mean(cr(hold)));
fprintf('Seasonal persistence MAE %.1f (model on same months %.1f)\n', mean(abs(naive(hm) - y(hm))), mean(abs(med(hm) - y(hm))));
tp = q(sum(ppc(hold,:), 1)); fprintf('Hold-out total: observed %.0f, predicted %.0f [%.0f, %.0f]\n', sum(y(hold)), tp(2), tp(1), tp(3));
% annual summaries
N = squeeze(sum(st(1:end-1,:,:), 2)); Se = squeeze(st(1:end-1,1,:)) + (ones(nm,1)*p.theta).*squeeze(st(1:end-1,2,:));
reff = (ones(nm,1)*p.R0).*Se./N; yrs = 2012:2026;
fprintf('\nYear  obs   pred-reported            infections(10^3)        Rc                 mean Reff\n');
RC = zeros(numel(idx), numel(yrs)); for d = 1:size(TH,2), RC(d,:) = est_rc_year(M, TH(:,d), yrs); end
ANN = zeros(numel(yrs), 14);
for k = 1:numel(yrs)
    j = M.D.year == yrs(k); a = q(sum(ppc(j & obs,:),1)); b = q(sum(inc(j,:),1)/1e3); c = q(RC(:,k)'); e = q(mean(reff(j,:),1));
    o = [2 1 3]; ANN(k,:) = [yrs(k) sum(y(j & obs)) a(o) b(o) c(o) e(o)];     % median, 2.5%, 97.5%
    fprintf('%d %5.0f  %5.0f [%5.0f,%5.0f]  %6.1f [%5.1f,%6.1f]  %.2f [%.2f,%.2f]  %.2f [%.2f,%.2f]\n', ANN(k,:));
end
% counterfactual SIAs
fprintf('\nInfections averted by each SIA (to May 2026), median [95%% CrI]:\n'); CF = zeros(numel(M.D.sia_frac), 3);
sub = TH(:, 1:min(NCF, size(TH,2))); inc_ref = inc(:, 1:size(sub,2));
for s = 1:numel(M.D.sia_frac)
    sf = M.D.sia_frac; sf(s) = 0; inc0 = est_simulate(M, sub, sf);
    t0 = M.D.datenum >= datenum(M.D.year(1),1,1) + (M.D.sia_idx(s) - 1);
    av = sum(inc0(t0,:) - inc_ref(t0,:), 1); CF(s,:) = q(av);
    fprintf('  SIA %d: %.1f [%.1f, %.1f] x 10^3\n', s, CF(s,2)/1e3, CF(s,1)/1e3, CF(s,3)/1e3);
end
C = corrcoef(flat);
save(fullfile(results_dir(), 'est_posterior.mat'), 'nat', 'ppc', 'lo', 'hi', 'lo50', 'hi50', 'med', 'inc', 'reff', 'ANN', 'CF', 'C', 'TH');
