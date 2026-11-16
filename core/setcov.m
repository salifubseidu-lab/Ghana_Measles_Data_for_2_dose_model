function Q = setcov(P, p1, p2, w1)
%SETCOV Set vaccination rates from uptake by the lifetime mapping:
%   nu1 = p1*mu/(1-p1), nu2 = p2*(ell+mu)/(1-p2); optional omega1 -> theta.
Q = P;
if nargin > 1 && ~isempty(p1), Q.nu1 = p1/(1 - p1)*P.mu; end
if nargin > 2 && ~isempty(p2), Q.nu2 = p2/(1 - p2)*(Q.ell + P.mu); end
if nargin > 3 && ~isempty(w1), Q.theta = 1 - w1; end
end
