function q = qhat_branch(lam, P)
%QHAT_BRANCH beta*mu on the equilibrium branch at force of infection lam:
%   the unique positive root of r q^2 + [r(c+d)lam - A - alpha] q + lam(r d c lam - A c - alpha d) = 0.
k = derived(P); th = P.theta;
A = k.A; al = P.alpha; d = P.d; c = d + P.K*P.Lambda;
G = (th*lam.^2 + (k.k2 + th*k.k1)*lam + k.D)./(k.Phi0 + th*lam);
r = P.eps./(k.k3*G);
a2 = r; a1 = r*(c + d).*lam - A - al; a0 = r*d*c.*lam.^2 - A*c*lam - al*d*lam;
disc = sqrt(a1.^2 - 4*a2.*a0);
q = zeros(size(lam));
pos = a1 <= 0;
q(pos)  = (-a1(pos) + disc(pos))./(2*a2(pos));
q(~pos) = (-2*a0(~pos))./(a1(~pos) + disc(~pos) + 1e-300);   % stable form
end
