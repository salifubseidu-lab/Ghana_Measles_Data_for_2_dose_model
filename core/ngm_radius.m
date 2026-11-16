function r = ngm_radius(P)
%NGM_RADIUS Spectral radius of F*V^{-1} computed numerically at the DFE.
x0 = dfe_state(P); J = model_jac(x0, P); N0 = sum(x0);
Fm = [0, P.beta*(x0(1) + P.theta*x0(2))/N0; 0 0];
Vm = -J([4 5],[4 5]) + Fm;
r = max(abs(eig(Fm/Vm)));
end
