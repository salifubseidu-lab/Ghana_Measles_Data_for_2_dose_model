function X = outbreak_ic(P, Rt0, I0, renorm)
%OUTBREAK_IC Endemic state with susceptibles raised so that the unsaturated
%   effective reproduction number equals Rt0, and I(0) = I0.  With renorm = true,
%   R is reset so that N equals the equilibrium population (used in Section 11).
if nargin < 2, Rt0 = 1.8; end
if nargin < 3, I0 = 10; end
if nargin < 4, renorm = false; end
X = upper_ee(P); N = sum(X); k = derived(P);
target = Rt0*N*k.k3*k.k40/(P.beta*P.eps);
dS = target - (X(1) + P.theta*X(2));  X(1) = X(1) + dS;  X(6) = X(6) - dS;
X(4) = I0*k.k40/P.eps;  X(5) = I0;
if renorm, X(6) = N - sum(X(1:5)); end
end
