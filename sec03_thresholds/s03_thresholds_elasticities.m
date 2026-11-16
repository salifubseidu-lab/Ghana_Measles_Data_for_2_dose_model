%S03_THRESHOLDS_ELASTICITIES  Derived quantities of Table 7 and elasticities of Table 9.
%   Reproduces: Rc, Q, Gamma, Kc, 1/Kc, nu1*, p1*, Rc(nu1->inf), required nu2,
%   R_SN for K/Kc = 2,5,10,100,150 and K->inf; elasticities Upsilon_x of Rc.
setup_paths();
P = baseline_params(); kc = Kc_fun(P); k = derived(P);
Rc = Rc_fun(P); Q = P.mu*k.Phi0/k.D;
fprintf('Rc = %.4f   Q = %.4f   Gamma = %.4f\n', Rc, Q, Gamma_fun(P));
fprintf('Kc = %.4e   1/Kc = %.1f per 10^6\n', kc, 1/kc);
% critical first-dose rate nu1* (Rc = 1) and p1*
f = @(nu1) Rc_fun(setfield(P, 'nu1', nu1)) - 1;            %#ok<SFLD>
nu1s = fzero(f, [P.nu1, 100]);
fprintf('nu1* = %.4f   p1* = %.4f\n', nu1s, nu1s/(nu1s + P.mu));
fprintf('Rc(nu1 -> inf) = R0*theta*mu/(mu+nu2) = %.4f\n', R0_fun(P)*P.theta*P.mu/(P.mu + P.nu2));
fprintf('required nu2 > mu*(theta*R0 - 1) = %.4f\n', P.mu*(P.theta*R0_fun(P) - 1));
% R_SN
krs = [2 5 10 100 150 1e6]; rsn = zeros(size(krs));
for i = 1:numel(krs), Q2 = P; Q2.K = krs(i)*kc; rsn(i) = RSN_fun(Q2); end
fprintf('R_SN for K/Kc = %s : %s\n', mat2str(krs), mat2str(rsn, 4));
fprintf('lower bound (gamma+d+mu)/k40 = %.4f\n', k.A/k.k40);
% elasticities Upsilon_x = (dRc/dx)(x/Rc), central differences
names = {'beta','gamma','nu1','alpha','nu2','theta','ell','d','eps','mu'};
ups = zeros(size(names)); hrel = 1e-6;
for i = 1:numel(names)
    nm = names{i}; Pp = P; Pm = P; x = P.(nm);
    Pp.(nm) = x*(1 + hrel); Pm.(nm) = x*(1 - hrel);
    if strcmp(nm, 'mu')        % Lambda = mu*N0 moves with mu
        Pp.Lambda = P.Lambda*(1 + hrel); Pm.Lambda = P.Lambda*(1 - hrel);
    end
    ups(i) = (Rc_fun(Pp) - Rc_fun(Pm))/(2*hrel*Rc);
end
for i = 1:numel(names), fprintf('  Upsilon_%-6s = %+.4f\n', names{i}, ups(i)); end
fid = fopen(fullfile(results_dir(), 'tab09_elasticities.csv'), 'w');
fprintf(fid, 'parameter,elasticity\n');
for i = 1:numel(names), fprintf(fid, '%s,%.6f\n', names{i}, ups(i)); end
fclose(fid);
