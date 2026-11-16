function out = prcc_outputs(P)
%PRCC_OUTPUTS [Rc, endemic incidence, endemic deaths, 10-year cases, 10-year deaths]
%   after 10 introductions at the DFE (cumulative cases/deaths integrated as extra states).
rc = Rc_fun(P); X = upper_ee(P);
if isempty(X), inc = 0; dth = 0; else, inc = P.eps*X(4); dth = P.d*X(5); end
X0 = dfe_state(P); X0(5) = 10; X0(4) = 10*(P.gamma + P.alpha + P.d + P.mu)/P.eps; X0(6) = 1e-3;
f = @(t,z) [model_rhs(z(1:6), P); P.eps*z(4); P.d*z(5)];
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'NonNegative', 1:8);
[~, Z] = ode_stiff(f, [0 5 10], [X0; 0; 0], opts);
out = [rc, inc, dth, Z(end,7), Z(end,8)];
end
