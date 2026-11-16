function P = random_params(mode)
%RANDOM_PARAMS Random parameter set ('measles' range or 'wide' range), as in Section A.
if strcmp(mode, 'measles')
    m = 1/(40 + 40*rand); N0 = 1e6;
    P.Lambda = m*N0; P.mu = m; P.nu1 = m*10^(-1 + 2.7*rand); P.ell = 0.2*rand;
    P.theta = 0.5*rand; P.eps = 365/(7 + 7*rand); P.gamma = 365/(5 + 5*rand);
    P.alpha = P.gamma*2*rand; p2 = 0.05 + 0.9*rand; P.nu2 = p2/(1 - p2)*(P.ell + m);
    cfr = 10^(-3 + (log10(0.15) + 3)*rand); P.d = cfr/(1 - cfr)*(P.gamma + P.alpha + m);
else
    m = 10^(-2 + 2*rand); N0 = 10^(2 + 4*rand);
    P.Lambda = m*N0; P.mu = m; P.nu1 = 10^(-2 + 4*rand); P.ell = 10^(-3 + 4*rand);
    P.theta = rand; P.eps = 10^(-1 + 3*rand); P.gamma = 10^(-1 + 3*rand);
    P.alpha = 10^(-1 + 3*rand); P.nu2 = 10^(-2 + 4*rand); P.d = 10^(-3 + 4.5*rand);
end
P.K = 0; P.beta = 1; kc = Kc_fun(P);
if rand < 0.85, P.K = max(kc, 1e-12)*10^(-1.5 + 3*rand); end
P.beta = beta_for_Rc(P, 10^(log10(0.3) + (log10(5) - log10(0.3))*rand));
end
