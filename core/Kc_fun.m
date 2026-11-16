function kc = Kc_fun(P)
%KC_FUN Treatment-capacity threshold Kc = k40*(k40*Gamma*k3/eps - d)/(alpha*Lambda).
k = derived(P);
kc = k.k40*(k.k40*Gamma_fun(P)*k.k3/P.eps - P.d)/(P.alpha*P.Lambda);
end
