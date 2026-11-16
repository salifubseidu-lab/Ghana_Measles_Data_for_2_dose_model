function r = Rc_fun(P)
%RC_FUN Control reproduction number (Theorem 3.3): Rc = R0*mu*Phi0/D.
k = derived(P); r = P.beta*P.mu*P.eps*k.Phi0/(k.k3*k.k40*k.D);
end
