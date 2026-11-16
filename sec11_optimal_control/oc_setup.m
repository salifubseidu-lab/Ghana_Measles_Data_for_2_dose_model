function S = oc_setup(P, X0, T, h, active, umax)
%OC_SETUP Discretisation and bookkeeping for the optimal-control problem.
S.P = P; S.X0 = X0(:); S.n = round(T/h); S.h = T/S.n; S.t = linspace(0, T, S.n + 1)';
S.active = active; S.umax = umax(:)'.*active; S.idx = find(active);
S.w = S.h*ones(S.n + 1, 1); S.w([1 end]) = S.h/2;   % trapezoid weights
S.scale = 1e-6;
end
