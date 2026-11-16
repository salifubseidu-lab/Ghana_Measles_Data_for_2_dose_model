function X = ee_from_lambda(lam, P)
%EE_FROM_LAMBDA Endemic equilibrium state for force of infection lam (Theorem 5.1).
k = derived(P); m = P.mu; th = P.theta;
I  = P.Lambda*lam/(P.beta*m + P.d*lam);
S  = P.Lambda*(k.k2 + th*lam)/((lam + k.k1)*(k.k2 + th*lam) - P.ell*P.nu1);
V1 = P.nu1*S/(k.k2 + th*lam);
E  = lam*(S + th*V1)/k.k3;
V2 = P.nu2*V1/m;
R  = (P.gamma + P.alpha/(1 + P.K*I))*I/m;
X  = [S; V1; V2; E; I; R];
end
