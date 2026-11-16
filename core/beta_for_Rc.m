function b = beta_for_Rc(P, r)
%BETA_FOR_RC Transmission rate giving control reproduction number r.
k = derived(P); b = r*k.k3*k.k40*k.D/(P.mu*P.eps*k.Phi0);
end
