%S10_BASIN  Finding 10.4 and Figure fig_n_basin: basin of the stable endemic equilibrium
%   in the bistable regime (K = 150 Kc, Rc = 0.9).  Initial states: susceptible fraction
%   sigma of the DFE level (rest immune) and I(0).  Grid sizes NA (full) and NB (zoom).
setup_paths();
if ~exist('NA', 'var'), NA = 14; end
if ~exist('NB', 'var'), NB = 36; end
P = baseline_params(); P.K = 150*Kc_fun(P); P.beta = beta_for_Rc(P, 0.9);
x0 = dfe_state(P); L = endemic_lambdas(P); Xu = ee_from_lambda(L(end), P); Xl = ee_from_lambda(L(1), P);
k = derived(P); sEE = (Xu(1) + Xu(2))/(x0(1) + x0(2)); sL = (Xl(1) + Xl(2))/(x0(1) + x0(2));
fprintf('stable EE: I* = %.1f, sigma = %.3f; saddle: I = %.3f, sigma = %.3f\n', Xu(5), sEE, Xl(5), sL);
grids = {linspace(0.04, 1, NA), logspace(-1, 4, NA); linspace(0.6, 0.95, NB), logspace(0, 4, NB)};
G = cell(1,2);
for g = 1:2
    sig = grids{g,1}; I0s = grids{g,2}; G{g} = zeros(numel(I0s), numel(sig));
    for i = 1:numel(I0s)
        for j = 1:numel(sig)
            X = x0; moved = (1 - sig(j))*(X(1) + X(2)); X(1) = X(1)*sig(j); X(2) = X(2)*sig(j);
            X(6) = moved + 1e-6; X(4) = I0s(i)*k.k40/P.eps; X(5) = I0s(i);
            G{g}(i,j) = basin_outcome(P, X);
        end
    end
    fprintf('grid %d: %d of %d initial states end endemic (%.1f%%)\n', g, sum(G{g}(:)), numel(G{g}), 100*mean(G{g}(:)));
end
[ii, jj] = find(G{2}); if ~isempty(ii)
    fprintf('zoom: endemic outcomes need sigma in [%.2f, %.2f], I(0) in [%.0f, %.0f]\n', min(grids{2,1}(jj)), max(grids{2,1}(jj)), min(grids{2,2}(ii)), max(grids{2,2}(ii)));
end
save(fullfile(results_dir(), 'basin.mat'), 'G', 'grids', 'sEE', 'sL');
fig = figure('Position', [100 100 950 370]);
for g = 1:2
    subplot(1,2,g); imagesc(grids{g,1}, log10(grids{g,2}), G{g}); set(gca, 'YDir', 'normal');
    colormap(flipud(gray)); caxis([0 1.3]); hold on; plot(sEE, log10(Xu(5)), 'rp', 'MarkerFaceColor', 'r');
    if g == 2, plot(sL, log10(Xl(5)), 'bo'); end
    xlabel('susceptible fraction \sigma of DFE level'); ylabel('log_{10} I(0)');
end
subplot(1,2,1); title('(a) full range (dark: endemic outcome)'); subplot(1,2,2); title('(b) zoom near the endemic state');
save_fig(fig, 'fig_n_basin');
