function h2 = h2_fun(lam, P)
%H2_FUN Analytic second derivative of h (Theorem 5.3, Step 2):
%   h'' = -2 pE q^2 A alpha K^2 Lambda^2/(k40 q + m lam)^3 - 2 theta^2 kappa/(theta lam + Phi0)^3,
%   q = beta*mu, c = d + K*Lambda, m = A c + alpha d, kappa = nu1 (nu2 + mu omega1).
k = derived(P); q = P.beta*P.mu; c = P.d + P.K*P.Lambda; m = k.A*c + P.alpha*P.d;
kap = P.nu1*(P.nu2 + P.mu*(1 - P.theta));
h2 = -2*k.pE*q^2*k.A*P.alpha*P.K^2*P.Lambda^2./(k.k40*q + m*lam).^3 ...
     - 2*P.theta^2*kap./(P.theta*lam + k.Phi0).^3;
end
