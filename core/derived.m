function k = derived(P)
%DERIVED Composite rates used throughout the paper.
k.k1  = P.nu1 + P.mu;
k.k2  = P.ell + P.nu2 + P.mu;
k.k3  = P.eps + P.mu;
k.k40 = P.gamma + P.alpha + P.d + P.mu;
k.D   = k.k1*k.k2 - P.ell*P.nu1;
k.Phi0 = k.k2 + P.theta*P.nu1;
k.A   = P.gamma + P.d + P.mu;
k.pE  = P.eps/k.k3;
k.pd  = P.d/k.k40;
end
