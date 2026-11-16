%S10_PRCC  Finding 10.5 and Figure fig_n_prcc: Latin hypercube sampling (NS samples) of
%   the ranges in Table 11 and partial rank correlation coefficients for five outputs.
setup_paths();
if ~exist('NS', 'var'), NS = 2000; end
rng(11);
names = {'R_0','latent','infectious','\alpha/\gamma','CFR','life exp.','p_1','p_2','\omega_1','\ell','log_{10}KN_0'};
lo = [12  8  6 0.0 0.01 50 0.60 0.30 0.75 0.00 3.0];
hi = [18 12 10 1.0 0.05 70 0.98 0.95 0.95 0.05 5.0];
X = lo + lhs_sample(NS, numel(lo)).*(hi - lo); Y = zeros(NS, 5);
for i = 1:NS
    Y(i,:) = prcc_outputs(prcc_params(X(i,:)));
    if mod(i, 200) == 0, fprintf('%d samples done\n', i); end
end
outs = {'R_c', 'endemic incidence', 'endemic deaths', '10-yr cases', '10-yr deaths'};
fprintf('%.1f%% of samples have Rc > 1\n', 100*mean(Y(:,1) > 1));
Rm = zeros(5, numel(lo)); Pm = Rm;
for o = 1:5
    [Rm(o,:), Pm(o,:)] = prcc_calc(X, Y(:,o));
    fprintf('%s:\n', outs{o});
    [~, ord] = sort(abs(Rm(o,:)), 'descend');
    for j = ord, fprintf('   %-14s %+.3f  (p = %.2g)\n', names{j}, Rm(o,j), Pm(o,j)); end
end
save(fullfile(results_dir(), 'prcc.mat'), 'X', 'Y', 'Rm', 'Pm', 'names', 'outs');
fig = figure('Position', [100 100 1100 380]);
for o = 1:5
    subplot(1,5,o); barh(Rm(o,:)); set(gca, 'YTick', 1:numel(names), 'YTickLabel', names, 'YDir', 'reverse');
    xlim([-1 1]); title(outs{o}); grid on;
end
save_fig(fig, 'fig_n_prcc');
