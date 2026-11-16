function v = dir_deriv(x, a, P, n)
%DIR_DERIV n-th directional derivative (n = 2 or 3) of the model vector field
%   at x along direction a (a may be complex):  d^n/ds^n F(x + s a) at s = 0.
%   Only the nonlinear terms contribute: T1 = beta*S*I/N, T2 = beta*theta*V1*I/N
%   and the treatment term g(I) = alpha*I/(1+K I).
N = sum(x); Na = sum(a);
T1 = rat_deriv(x(1), a(1), x(5), a(5), N, Na, n)*P.beta;
T2 = rat_deriv(x(2), a(2), x(5), a(5), N, Na, n)*P.beta*P.theta;
I = x(5);
gn = P.alpha*(-1)^(n+1)*factorial(n)*P.K^(n-1)/(1 + P.K*I)^(n+1);   % g^(n)(I)
G = gn*a(5)^n;
v = [-T1; -T2; 0; T1 + T2; -G; G];
end

function dn = rat_deriv(u, ua, w, wa, N, Na, n)
% n-th derivative at s=0 of (u+s ua)(w+s wa)/(N+s Na)
p = [u*w, u*wa + ua*w, ua*wa];
r = -Na/N;  an = 0;
for kk = 0:min(n,2)
    an = an + p(kk+1)*r^(n-kk)/N;
end
dn = factorial(n)*an;
end
