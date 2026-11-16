%S06_BIF2P_CODIM2  Two-parameter bifurcation diagram in the (K/Kc, Rc) plane
%   (Figure fig_bif2p) and the codimension-two points (Finding 6.11):
%   fold (saddle-node) curve, Hopf curves H1/H2 with the sign of l1,
%   Bogdanov-Takens point (bisection) and the two generalised-Hopf (Bautin) points.
setup_paths();
if ~exist('NK', 'var'), NK = 241; end           % number of K values on the coarse grid
P0 = baseline_params(); kc = Kc_fun(P0);
Ks = kc*logspace(-3, 3, NK); FOLD = zeros(0,3); HOPF = zeros(0,6);
for i = 1:numel(Ks)
    P = P0; P.K = Ks(i);
    [lsn, rsn] = fold_point(P);
    if ~isnan(lsn), FOLD(end+1,:) = [P.K/kc, rsn, lsn]; end %#ok<SAGROW>
    HOPF = [HOPF; hopf_points(P, 12, 500)]; %#ok<AGROW>
end
% extra resolution near (Kc, 1) for the tail of H1
for kr = [logspace(-1, 0, 25), linspace(1.0, 1.0139, 6)]
    P = P0; P.K = kr*kc; HOPF = [HOPF; hopf_points(P, 1.5, 900)]; %#ok<AGROW>
end
fprintf('fold points %d, Hopf points %d (l1>0: %d, l1<0: %d)\n', size(FOLD,1), size(HOPF,1), sum(HOPF(:,4)>0), sum(HOPF(:,4)<0));
% --- generalised Hopf points: zeros of l1 along H1 (first Hopf point at given K)
l1first = @(kr) first_l1(P0, kr*kc);
GH = zeros(2,6); br = [0.1413 0.1778; 0.6 0.75];
for g = 1:2
    kr = fzero(l1first, br(g,:), optimset('TolX', 1e-6));
    P = P0; P.K = kr*kc; H = hopf_points(P, 1.5, 800); GH(g,:) = H(1,:);
    fprintf('GH%d: K/Kc = %.4f, Rc = %.4f, omega = %.3f, I* = %.3f\n', g, kr, H(1,2), H(1,5), H(1,6));
end
% --- Bogdanov-Takens point: bisection on the existence of H1 above the fold
lo = 1.01; hi = 1.05;
for it = 1:40
    mid = (lo + hi)/2; P = P0; P.K = mid*kc;
    if ~isempty(hopf_points(P, 1.5, 1500)), lo = mid; else, hi = mid; end
end
P = P0; P.K = hi*kc; [lsn, rsn] = fold_point(P); Q = P_on_branch(lsn, P); X = ee_from_lambda(lsn, Q);
ev = sort(eig(model_jac(X, Q)));
BT = [hi, rsn, lsn, X(5)];
fprintf('BT: K/Kc = %.4f, R_SN = %.5f, I* = %.3f\n', BT(1), BT(2), BT(4));
fprintf('eigenvalues at the fold equilibrium (BT):'); fprintf(' %.3g%+.3gi', [real(ev) imag(ev)]'); fprintf('\n');
save(fullfile(results_dir(), 'bif2p.mat'), 'FOLD', 'HOPF', 'GH', 'BT');
% --- figure
fig = figure('Position', [100 100 950 380]);
for pan = 1:2
    subplot(1,2,pan); hold on;
    plot(FOLD(:,1), FOLD(:,2), 'k-', 'LineWidth', 1.2);
    s = HOPF(:,4) < 0;
    plot(HOPF(s,1), HOPF(s,2), '.', 'Color', [0.165 0.435 0.592], 'MarkerSize', 9);
    plot(HOPF(~s,1), HOPF(~s,2), '.', 'Color', [0.769 0.271 0.212], 'MarkerSize', 9);
    plot([1e-3 1e3], [1 1], '--', 'Color', [0.4 0.4 0.4]);
    plot(1, 1, 'ks'); plot(BT(1), BT(2), '^', 'Color', [0.737 0.424 0.145], 'MarkerFaceColor', [0.737 0.424 0.145]);
    plot(GH(:,1), GH(:,2), 'd', 'Color', [0.482 0.173 0.749], 'MarkerFaceColor', [0.482 0.173 0.749]);
    set(gca, 'XScale', 'log'); xlabel('K/K_c'); ylabel('R_c');
    if pan == 1
        set(gca, 'YScale', 'log'); xlim([5e-3 1e3]); ylim([0.65 12]); title('(a) (K/K_c, R_c) plane');
        legend({'SN (fold)', 'Hopf, l_1<0', 'Hopf, l_1>0', 'TC (R_c=1)', 'TC-SN', 'BT', 'GH (Bautin)'}, 'Location', 'north');
    else
        xlim([0.08 2]); ylim([0.995 1.03]); title('(b) zoom near (K_c, 1)');
    end
end
save_fig(fig, 'fig_bif2p');
