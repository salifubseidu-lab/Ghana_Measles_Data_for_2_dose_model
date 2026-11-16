function l = lam_for_rc(P, rc)
%LAM_FOR_RC Force of infection of the upper endemic equilibrium at Rc = rc.
Q = P; Q.beta = beta_for_Rc(P, rc); L = endemic_lambdas(Q); l = L(end);
end
