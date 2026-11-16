function x0 = dfe_state(P)
%DFE_STATE Disease-free equilibrium X0 (Lemma 3.2).
k = derived(P);
x0 = [P.Lambda*k.k2/k.D; P.Lambda*P.nu1/k.D; P.Lambda*P.nu1*P.nu2/(P.mu*k.D); 0; 0; 0];
end
