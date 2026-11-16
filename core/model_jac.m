function J = model_jac(x, P)
%MODEL_JAC Analytic Jacobian of model (2.1) at state x (6-by-1).
S = x(1); V1 = x(2); I = x(5);
N = sum(x);
b = P.beta; th = P.theta;
% transmission terms T1 = b*S*I/N, T2 = b*th*V1*I/N
e = ones(1,6); e5 = [0 0 0 0 1 0];
gT1 = b*( I/N*[1 0 0 0 0 0] + S/N*e5 - S*I/N^2*e );
gT2 = b*th*( I/N*[0 1 0 0 0 0] + V1/N*e5 - V1*I/N^2*e );
% treatment g(I) = alpha*I/(1+K I), g'(I) = alpha/(1+K I)^2
gp = P.alpha/(1 + P.K*I)^2;
J = zeros(6);
J(1,:) = -gT1;  J(1,1) = J(1,1) - (P.nu1 + P.mu);  J(1,2) = J(1,2) + P.ell;
J(2,:) = -gT2;  J(2,1) = J(2,1) + P.nu1;  J(2,2) = J(2,2) - (P.ell + P.nu2 + P.mu);
J(3,2) = P.nu2;  J(3,3) = -P.mu;
J(4,:) = gT1 + gT2;  J(4,4) = J(4,4) - (P.eps + P.mu);
J(5,4) = P.eps;  J(5,5) = -(P.gamma + gp + P.d + P.mu);
J(6,5) = P.gamma + gp;  J(6,6) = -P.mu;
end
