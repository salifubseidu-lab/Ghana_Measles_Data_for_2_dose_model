%S06_REGIME_MAP  Figure fig_regime: regimes in the (K/Kc, Rc) plane (local analysis of
%   the upper equilibrium), and R_SN by bisection for several K/Kc.
%   codes: 0 no EE; 1 two EE, upper stable; 2 two EE, upper unstable;
%          3 unique EE stable; 4 unique EE unstable (cycles)
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0);
Kr = logspace(-3, 2.5, 221); Rr = linspace(0.3, 6, 229); code = zeros(numel(Rr), numel(Kr));
for j = 1:numel(Kr)
    for i = 1:numel(Rr)
        P = P0; P.K = Kr(j)*kc; P.beta = beta_for_Rc(P, Rr(i)); [X, ev, L] = upper_ee(P);
        if isempty(L), code(i,j) = 0; continue; end
        st = max(real(ev)) < 0;
        if numel(L) == 2, code(i,j) = 1 + ~st; else, code(i,j) = 3 + ~st; end
    end
end
for c = 0:4, fprintf('code %d: %d cells\n', c, sum(code(:) == c)); end
for kr = [1.5 2 5 10 40 100]
    P = P0; P.K = kr*kc; fprintf('K/Kc = %g: R_SN = %.5f\n', kr, RSN_fun(P));
end
save(fullfile(results_dir(), 'regime.mat'), 'Kr', 'Rr', 'code');
cmap = [0.949 0.949 0.949; 0.659 0.855 0.863; 0.898 0.596 0.608; 0.271 0.482 0.616; 0.914 0.769 0.416];
fig = figure('Position', [100 100 640 420]);
imagesc(log10(Kr), Rr, code); set(gca, 'YDir', 'normal'); colormap(cmap); caxis([-0.5 4.5]);
hold on; plot([0 0], [0.3 6], 'k:'); plot(log10(Kr([1 end])), [1 1], 'k:');
xlabel('log_{10}(K/K_c)'); ylabel('R_c');
cb = colorbar; labs = {'no EE','DFE + stable EE','DFE + unstable EE','unique stable EE','unique unstable EE'};
try, set(cb, 'Ticks', 0:4, 'TickLabels', labs);          % MATLAB R2014b+
catch, set(cb, 'YTick', 0:4, 'YTickLabel', labs); end    % Octave
save_fig(fig, 'fig_regime');
