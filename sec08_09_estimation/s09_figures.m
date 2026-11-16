%S09_FIGURES  Figures fig_e_fit, fig_e_post, fig_e_prof, fig_e_corr/trace and fig_e_sens.
setup_paths(); rd = results_dir();
M = est_setup('covid', true); P = load(fullfile(rd, 'est_posterior.mat')); S = load(fullfile(rd, 'est_mcmc.mat'));
dates = M.D.datenum; y = M.D.y; fit = M.D.fit; C1 = [0.165 0.435 0.592]; C2 = [0.769 0.271 0.212]; C3 = [0.416 0.6 0.306];
p = est_unpack(M, P.TH); meanrep = prctile_((P.inc.*(ones(size(P.inc,1),1)*p.rho))', 50)';
fig = figure('Position', [100 100 950 620]);
subplot(3,1,[1 2]); hold on;
fill([dates; flipud(dates)], [P.lo; flipud(P.hi)], C1, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
fill([dates; flipud(dates)], [P.lo50; flipud(P.hi50)], C1, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
plot(dates, meanrep, 'Color', C1, 'LineWidth', 1.2);
plot(dates(fit), y(fit), 'k.', 'MarkerSize', 9); plot(dates(~fit), y(~fit), 'o', 'Color', C2, 'MarkerSize', 4);
sd = M.D.days(M.D.sia_idx); for s = 1:numel(sd), plot([sd(s) sd(s)], [0 3000], '--', 'Color', [0.4 0.4 0.4]); end
ylim([0 3000]); ylabel('cases per month'); datetick('x', 'yyyy'); xlim([dates(1) dates(end)]); legend({'95% predictive', '50% predictive', 'median \rho C_m', 'observed (fit)', 'observed (hold-out)'}, 'Location', 'northwest');
subplot(3,1,3); hold on; r = prctile_(P.reff', [2.5 50 97.5])'; d2 = dates(1:size(r,1));
fill([d2; flipud(d2)], [r(:,1); flipud(r(:,3))], C3, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); plot(d2, r(:,2), 'Color', C3);
plot(d2([1 end]), [1 1], 'k-'); ylabel('R_{eff} (seasonal mean)'); datetick('x', 'yyyy'); xlim([dates(1) dates(end)]);
save_fig(fig, 'fig_e_fit');
% posteriors vs priors
fl = S.flat; rho = 1./(1 + exp(-fl(:,5)));
pan = {fl(:,1), 'R_0'; fl(:,7), '\omega_1'; log10(rho), 'log_{10}\rho'; log10(fl(:,9)), 'log_{10}\lambda_{ref}'; fl(:,6), 's_0';
       fl(:,2), 'A_s'; fl(:,4), 'log_{10}\iota'; fl(:,10), 'c_{npi}'; fl(:,11), 'h_{SIA}'; exp(fl(:,8)), '\phi_{NB}'};
fig = figure('Position', [100 100 1000 400]);
for i = 1:size(pan,1)
    subplot(2,5,i); hist(pan{i,1}, 50); xlabel(pan{i,2}); set(gca, 'YTick', []);
end
save_fig(fig, 'fig_e_post');
% profiles
if exist(fullfile(rd, 'est_profiles.mat'), 'file')
    F = load(fullfile(rd, 'est_profiles.mat')); fn = fieldnames(F.PROF); fig = figure('Position', [100 100 1100 260]);
    for i = 1:numel(fn)
        subplot(1, numel(fn), i); g = F.PROF.(fn{i}).grid; if strcmp(fn{i}, 'logitrho'), g = 1./(1 + exp(-g)); end
        plot(g, F.PROF.(fn{i}).delta, 'o-', 'Color', C1); hold on; plot(g([1 end]), [1.92 1.92], '--', 'Color', C2);
        ylim([0 8]); xlabel(fn{i});
        if any(strcmp(fn{i}, {'kappa', 'logitrho'})), set(gca, 'XScale', 'log'); end
    end
    save_fig(fig, 'fig_e_prof');
end
% correlation and trace
fig = figure('Position', [100 100 900 360]);
subplot(1,2,1); imagesc(P.C); caxis([-1 1]); colorbar; axis square; set(gca, 'XTick', 1:M.npar, 'YTick', 1:M.npar, 'XTickLabel', M.names, 'YTickLabel', M.names);
title('posterior correlation');
if ~isempty(S.lp), subplot(1,2,2); plot(S.lp(1:20:end,:), 'Color', [0.70 0.80 0.88]); xlabel('iteration / 20'); ylabel('log posterior'); end
save_fig(fig, 'fig_e_corr');
% sensitivity
if exist(fullfile(rd, 'est_sensitivity.mat'), 'file')
    E = load(fullfile(rd, 'est_sensitivity.mat')); fig = figure('Position', [100 100 1000 380]); cols = [1 2 3 7]; lab = {'R_0', '\omega_1', '\rho', 'R_c (2024)'};
    for i = 1:4
        subplot(1,4,i); plot(E.T(:,cols(i)), size(E.T,1):-1:1, 'o', 'Color', C1); hold on;
        plot(E.T(1,cols(i))*[1 1], [0 size(E.T,1)+1], '--', 'Color', C2); xlabel(lab{i});
        if i == 1, set(gca, 'YTick', 1:size(E.T,1), 'YTickLabel', flipud(E.sc(:,1))); else, set(gca, 'YTick', []); end
    end
    save_fig(fig, 'fig_e_sens');
end
