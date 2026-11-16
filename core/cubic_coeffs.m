function c = cubic_coeffs(P)
%CUBIC_COEFFS Coefficients [c3 c2 c1 c0] of the cubic P(lambda) of Theorem 5.1
%   P(lam) = k3*[k40*beta*mu + M*lam]*Psi(lam)
%          - eps*(beta*mu + d*lam)*(beta*mu + (d+K*Lambda)*lam)*Phi(lam).
k = derived(P); m = P.mu; th = P.theta; bm = P.beta*m;
M  = (P.gamma + P.d + m)*(P.d + P.K*P.Lambda) + P.alpha*P.d;
A1 = [M, k.k40*bm];
Psi = [th, k.k2 + th*k.k1, k.D];
B1 = [P.d, bm];  Ph = [th, k.Phi0];  C1 = [P.d + P.K*P.Lambda, bm];
c = k.k3*conv(A1, Psi) - P.eps*conv(conv(B1, Ph), C1);
end
