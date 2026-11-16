function X = oc_forward(U, X0, P, h)
%OC_FORWARD RK4 state sweep; controls at the half step are averaged from the grid.
n = size(U,1) - 1; X = zeros(n+1, 6); X(1,:) = X0(:)';
for k = 1:n
    x = X(k,:)'; um = 0.5*(U(k,:) + U(k+1,:))';
    k1 = oc_f(x, U(k,:)', P); k2 = oc_f(x + h/2*k1, um, P); k3 = oc_f(x + h/2*k2, um, P); k4 = oc_f(x + h*k3, U(k+1,:)', P);
    X(k+1,:) = (x + h/6*(k1 + 2*k2 + 2*k3 + k4))';
end
end
