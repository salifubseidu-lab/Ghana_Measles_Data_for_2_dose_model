function f = oc_f(x, u, P)
%OC_F Controlled model (11.1); x = [S V1 V2 E I R]', u = [u1 u2 u3]'.
S = x(1); V1 = x(2); V2 = x(3); E = x(4); I = x(5); R = x(6); N = S + V1 + V2 + E + I + R;
lam = P.beta*I/N; tau = P.alpha*u(3) + P.alpha*(1 - u(3))/(1 + P.K*I);
f = [P.Lambda + P.ell*V1 - lam*S - (P.nu1 + u(1) + P.mu)*S;
     (P.nu1 + u(1))*S - P.theta*lam*V1 - (P.ell + P.nu2 + u(2) + P.mu)*V1;
     (P.nu2 + u(2))*V1 - P.mu*V2;
     lam*(S + P.theta*V1) - (P.eps + P.mu)*E;
     P.eps*E - (P.gamma + tau + P.d + P.mu)*I;
     (P.gamma + tau)*I - P.mu*R];
end
