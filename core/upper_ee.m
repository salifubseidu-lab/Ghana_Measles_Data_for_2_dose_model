function [X, ev, lams] = upper_ee(P)
%UPPER_EE Upper endemic equilibrium and its Jacobian eigenvalues ([] if none).
lams = endemic_lambdas(P); X = []; ev = [];
if isempty(lams), return; end
X = ee_from_lambda(lams(end), P); ev = eig(model_jac(X, P));
end
