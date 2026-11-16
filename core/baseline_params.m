function P = baseline_params()
%BASELINE_PARAMS Illustrative baseline of Table 7 (time unit: years).
%   P is a struct with fields Lambda, mu, nu1, nu2, ell, theta, eps, gamma,
%   alpha, d, K, beta.  beta is chosen so that R0 = 15, K = 0.
mu = 1/60; N0 = 1e6;
P.Lambda = mu*N0;  P.mu = mu;
P.nu1 = 4*mu;                     % p1 = 0.8
P.ell = 0.01;  P.theta = 0.15;    % omega1 = 0.85
P.nu2 = 1.5*(P.ell + mu);         % p2 = 0.6
P.eps = 365/10;  P.gamma = 365/10;
P.alpha = 365/7 - 365/10;         % 7-day infectious period when unsaturated
P.d = 0.03/0.97*(P.gamma + P.alpha + mu);   % CFR 3% under unsaturated treatment
P.K = 0;
P.beta = 1;
k = derived(P);
P.beta = 15*k.k3*k.k40/P.eps;     % R0 = 15
end
