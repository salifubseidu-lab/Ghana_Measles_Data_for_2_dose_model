function [J, g, X, Z] = oc_objective(zv, S)
%OC_OBJECTIVE Objective J and its adjoint gradient with respect to the scaled controls
%   zv in [0,1] (only active controls are optimisation variables).  S is the setup
%   struct from OC_SETUP.  The returned J and g are scaled by S.scale.
U = oc_unpack(zv, S);
X = oc_forward(U, S.X0, S.P, S.h); Z = oc_backward(X, U, S.P, S.h);
J = trapz(S.t, oc_lagr(X, U, S.P));
G = oc_dHdu(X, U, Z, S.P);
g = G(:, S.idx).*(S.w*S.umax(S.idx)); g = g(:);
J = J*S.scale; g = g*S.scale;
end
