%S11_COST_EFFECTIVENESS  Table 16 and Figure fig_cea: ACER, dominance and ICERs per case
%   and per death averted (intervention cost only), from the late-detection results.
setup_paths();
S = load(fullfile(results_dir(), 'oc_results.mat')); L = S.RES.late;
keys = {'s000','s100','s010','s001','s110','s101','s011','s111'}; nm = {'none','A','B','C','D','E','F','G'};
fig = figure('Position', [100 100 900 360]);
for e = 1:2
    fld = {'cases','deaths'}; f = fld{e};
    eff = cellfun(@(k) L.s000.(f) - L.(k).(f), keys); cost = cellfun(@(k) L.(k).cost, keys);
    T = cea_frontier(nm, eff, cost);
    fprintf('\nper %s averted:\n', f(1:end-1));
    for i = 1:numel(T)
        st = T(i).status; if isempty(st), st = 'frontier'; end
        fprintf('  %-5s averted %9.0f  cost %9.0f  ACER %8.2f  ICER %8.2f  %s\n', T(i).key, T(i).eff, T(i).cost, T(i).acer, T(i).icer, st);
    end
    subplot(1,2,e); hold on; fr = strcmp({T.status}, '');
    plot([T(fr).eff], [T(fr).cost], 'k-o', 'MarkerFaceColor', 'k'); plot([T(~fr).eff], [T(~fr).cost], 'ko');
    for i = 1:numel(T), text(T(i).eff, T(i).cost, ['  ' T(i).key]); end
    xlabel(sprintf('%s averted', f)); ylabel('intervention cost (US$)');
end
net = cellfun(@(k) L.(k).cost - 150*(L.s000.cases - L.(k).cases), keys(2:end));
fprintf('\nnet cost (cost - w x cases averted): %s\n', mat2str(round(net)));
save_fig(fig, 'fig_cea');
