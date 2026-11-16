function g = Gamma_fun(P)
%GAMMA_FUN Gamma = 1 - theta*nu1*(nu2 + mu*omega1)/Phi0^2  (Lemma 5.2).
k = derived(P);
g = 1 - P.theta*P.nu1*(P.nu2 + P.mu*(1 - P.theta))/k.Phi0^2;
end
