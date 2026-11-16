function g = cx_maxre(lam, P)
%CX_MAXRE Largest real part among complex eigenvalues at the branch point lam
%   (NaN if none); used to locate Hopf points along the branch.
Q = P_on_branch(lam, P); X = ee_from_lambda(lam, Q); ev = eig(model_jac(X, Q));
cx = ev(abs(imag(ev)) > 1e-12);
if isempty(cx), g = NaN; else, g = max(real(cx)); end
end
