%S09_PROFILES  Practical identifiability (Section 9.6): profile posteriors of R0,
%   lambda_ref, omega1, rho (logit), s0 and h_SIA by continuation (Nelder-Mead from the
%   MAP and from the neighbouring optimum), with the profiled parameter's prior removed;
%   95% interval where the profile is within 1.92 of its minimum (Raue et al. 2009).
setup_paths();
if ~exist('PMAXEVAL', 'var'), PMAXEVAL = 4000; end   % evaluations per Nelder-Mead run
S = load(fullfile(results_dir(), 'est_mcmc.mat')); M = est_setup('covid', true); xmap = S.R.x;
if ~exist('PSPECS', 'var')
specs = {1, [6 7 8 9 10 11 12 14 16 18]; 9, [0.0002 0.001 0.002 0.003 0.004 0.006 0.008 0.015 0.03 0.06 0.12 0.2 0.3];
         7, [0.55 0.65 0.72 0.78 0.84 0.9 0.95]; 5, [-5 -4 -3 -2.2 -1.5 -0.8 0 0.8 1.6];
         6, [0.01 0.02 0.03 0.04 0.05 0.07 0.1]; 11, [0.1 0.25 0.4 0.55 0.7 0.85 1.0 1.2 1.5]};
else, specs = PSPECS; end
PROF = struct();
for s = 1:size(specs,1)
    j = specs{s,1}; grid = specs{s,2}; nm = M.names{j};
    obj = @(v, val) -profile_obj(M, j, val, v);
    x0 = xmap; x0(j) = []; vals = zeros(size(grid)); X = zeros(numel(xmap), numel(grid));
    up = find(grid >= xmap(j)); dn = fliplr(find(grid < xmap(j)));
    for seq = {up, dn}
        prev = x0;
        for g = seq{1}
            f = @(v) obj(v(:), grid(g)); best = Inf;
            for st = {prev, x0}
                [v, fv] = nm_polish(f, st{1}, PMAXEVAL, 6, 1e-7);   % Nelder-Mead with restarts
                if fv < best, best = fv; vb = v(:); end
            end
            prev = vb; vals(g) = best; X(:,g) = [vb(1:j-1); grid(g); vb(j:end)];
            fprintf('%s = %g: %.3f\n', nm, grid(g), best);
        end
    end
    fmap = -profile_obj(M, j, xmap(j), x0);               % value at the MAP itself (own prior removed)
    d = vals - min([vals, fmap]); gg = grid;
    if j == 5, gg = 1./(1 + exp(-grid)); end      % interpolate reporting on the probability scale
    ci = profile_ci(gg, d);
    fprintf('  95%% profile interval for %s: [%s, %s]\n', nm, num2str(ci(1), 4), num2str(ci(2), 4));
    PROF.(nm) = struct('grid', grid, 'delta', d, 'X', X, 'ci', ci);
end
save(fullfile(results_dir(), 'est_profiles.mat'), 'PROF');
