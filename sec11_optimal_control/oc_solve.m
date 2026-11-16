function R = oc_solve(P, X0, active, varargin)
%OC_SOLVE Solve the optimal-control problem for the active controls (logical 1x3),
%   from two starting guesses (u = 0 and u = umax/2), keeping the better solution.
%   Options (name/value): 'T' (5), 'h' (0.002), 'umax' ([12 6 0.9]), 'maxit' (3000).
o = struct('T', 5, 'h', 0.002, 'umax', [12 6 0.9], 'maxit', 3000, 'starts', [0 0.5]);
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
S = oc_setup(P, X0, o.T, o.h, active, o.umax); nv = (S.n + 1)*numel(S.idx);
best = []; 
if isempty(S.idx), o.starts = 0; end
for st = o.starts
    if isempty(S.idx), zv = zeros(0,1); info = struct('it', 0);
    else, [zv, info] = oc_pqn(@(z) oc_objective(z, S), st*ones(nv,1), o.maxit); end
    U = oc_unpack(zv, S); X = oc_forward(U, S.X0, P, S.h); Z = oc_backward(X, U, P, S.h);
    J = trapz(S.t, oc_lagr(X, U, P));
    if isempty(best) || J < best.J, best = struct('J', J, 'U', U, 'X', X, 'Z', Z, 'start', st, 'iters', info.it); end
end
R = best; R.t = S.t; R.S = S;
% outcomes
X = R.X; U = R.U; t = R.t; I = X(:,5); extra = U(:,3).*P.alpha*P.K.*I.^2./(1 + P.K*I);
R.cases = trapz(t, P.eps*X(:,4)); R.deaths = trapz(t, P.d*I);
R.cost = trapz(t, P.c1*U(:,1).*X(:,1) + P.c2*U(:,2).*X(:,2) + P.c3*extra + 0.5*(P.b1*U(:,1).^2 + P.b2*U(:,2).^2 + P.b3*U(:,3).^2));
R.peakI = max(I); R.doses1 = trapz(t, U(:,1).*X(:,1)); R.doses2 = trapz(t, U(:,2).*X(:,2)); R.courses = trapz(t, extra);
Us = oc_ustar(X, R.Z, P, S.umax); a = logical(active);
if any(a), R.fp = max(max(abs(U(:,a) - Us(:,a))./(ones(size(U,1),1)*S.umax(a)))); else, R.fp = 0; end
end
