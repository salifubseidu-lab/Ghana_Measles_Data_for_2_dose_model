function g = oc_adj(x, u, z, P)
%OC_ADJ Adjoint right-hand side  dzeta/dt = -dH/dX  (Theorem 11.4).
S = x(1); V1 = x(2); I = x(5); N = sum(x); lam = P.beta*I/N; th = P.theta;
zS = z(1); zV1 = z(2); zV2 = z(3); zE = z(4); zI = z(5); zR = z(6);
Q = P.beta*I/N^2*(S*(zE - zS) + th*V1*(zE - zV1));
KI = P.K*I;
dS  = P.c1*u(1) + (P.nu1 + u(1))*(zV1 - zS) - P.mu*zS + lam*(zE - zS) - Q;
dV1 = P.c2*u(2) + P.ell*(zS - zV1) + (P.nu2 + u(2))*(zV2 - zV1) - P.mu*zV1 + th*lam*(zE - zV1) - Q;
dV2 = -P.mu*zV2 - Q;
dE  = P.w*P.eps + P.eps*(zI - zE) - P.mu*zE - Q;
dI  = P.c3*u(3)*P.alpha*KI*(2 + KI)/(1 + KI)^2 + P.beta/N*(S*(zE - zS) + th*V1*(zE - zV1)) - Q ...
      + (P.gamma + P.alpha*u(3) + P.alpha*(1 - u(3))/(1 + KI)^2)*(zR - zI) - (P.d + P.mu)*zI;
dR  = -P.mu*zR - Q;
g = -[dS; dV1; dV2; dE; dI; dR];
end
