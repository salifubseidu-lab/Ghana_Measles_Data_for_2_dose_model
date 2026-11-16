function P = prcc_params(x)
%PRCC_PARAMS Parameter set from an LHS point x = [R0, latent (d), infectious (d),
%   alpha/gamma, CFR, life expectancy (yr), p1, p2, omega1, ell, log10(K N0)] (Table 11).
mu = 1/x(6); N0 = 1e6;
P.Lambda = mu*N0; P.mu = mu; P.eps = 365/x(2); P.gamma = 365/x(3); P.ell = x(10); P.theta = 1 - x(9);
P.alpha = x(4)*P.gamma; P.nu1 = x(7)/(1 - x(7))*mu; P.nu2 = x(8)/(1 - x(8))*(x(10) + mu);
P.d = x(5)/(1 - x(5))*(P.gamma + P.alpha + mu); P.K = 10^x(11)/N0; P.beta = 1;
k = derived(P); P.beta = x(1)*k.k3*k.k40/P.eps;
end
