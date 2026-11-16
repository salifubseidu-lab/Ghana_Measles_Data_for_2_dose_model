function h = h_fun(lam, P)
%H_FUN The scalar function h(lambda) = 1/L(lambda) - Psi(lambda)/Phi(lambda) (Theorem 5.1).
k = derived(P); th = P.theta; bm = P.beta*P.mu;
c = P.d + P.K*P.Lambda;
L = (k.k3/P.eps)*(k.A./(bm + P.d*lam) + P.alpha./(bm + c*lam));
h = 1./L - (th*lam.^2 + (k.k2 + th*k.k1)*lam + k.D)./(k.Phi0 + th*lam);
end
