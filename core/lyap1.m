function [l1, om, re, q] = lyap1(x, P)
%LYAP1 First Lyapunov coefficient (Kuznetsov 2004, eq. 3.20) at equilibrium x.
%   Uses the complex pair with smallest |Re|.  Returns l1, omega, Re of the pair and
%   the unit right eigenvector q (needed for the normal-form amplitude).
A = model_jac(x, P);
[V, Dg] = eig(A); ev = diag(Dg);
score = abs(real(ev)) + 1e9*(abs(imag(ev)) < 1e-12) + 1e9*(imag(ev) < 0);
[~, idx] = min(score);
om = imag(ev(idx)); re = real(ev(idx)); q = V(:,idx);
[W, Dl] = eig(A.'); evl = diag(Dl);
[~, j] = min(abs(evl - (-1i*om)));
p = W(:,j);
p = p/conj(p'*q);                       % normalise so that <p,q> = p'*q = 1
qb = conj(q); I6 = eye(6);
t1 = p'*Cform(x, q, q, qb, P);
t2 = -2*(p'*Bform(x, q, A\Bform(x, q, qb, P), P));
t3 = p'*Bform(x, qb, (2i*om*I6 - A)\Bform(x, q, q, P), P);
l1 = real(t1 + t2 + t3)/(2*om);
end
