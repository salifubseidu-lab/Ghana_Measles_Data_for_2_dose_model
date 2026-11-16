function Jy = jac_log(~, y, P)
%JAC_LOG Jacobian of ODE_LOG.
x = exp(max(y, -700)); J = model_jac(x, P); f = model_rhs(x, P);
Jy = J.*(ones(6,1)*x.')./(x*ones(1,6)) - diag(f./x);
end
