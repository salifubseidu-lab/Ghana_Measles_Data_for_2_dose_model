function r = rho_branch(lam, P)
%RHO_BRANCH Branch function Rc = rho(lambda) of Theorem 6.3.
k = derived(P);
r = P.eps*k.Phi0*qhat_branch(lam, P)/(k.k3*k.k40*k.D);
end
