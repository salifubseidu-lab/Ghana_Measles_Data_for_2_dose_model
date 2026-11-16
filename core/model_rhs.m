function f = model_rhs(x, P)
%MODEL_RHS Right-hand side of model (2.1), x = [S;V1;V2;E;I;R] (column).
%   Works column-wise if x is 6-by-n.
S = x(1,:); V1 = x(2,:); V2 = x(3,:); E = x(4,:); I = x(5,:); R = x(6,:);
N = S + V1 + V2 + E + I + R;
lam = P.beta*I./N;
tr  = P.gamma + P.alpha./(1 + P.K*I);           % removal by recovery + saturated treatment
f = [P.Lambda + P.ell*V1 - lam.*S - (P.nu1 + P.mu)*S;
     P.nu1*S - P.theta*lam.*V1 - (P.ell + P.nu2 + P.mu)*V1;
     P.nu2*V1 - P.mu*V2;
     lam.*(S + P.theta*V1) - (P.eps + P.mu)*E;
     P.eps*E - (tr + P.d + P.mu).*I;
     tr.*I - P.mu*R];
end
