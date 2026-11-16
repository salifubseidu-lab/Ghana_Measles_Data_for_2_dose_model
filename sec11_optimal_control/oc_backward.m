function Z = oc_backward(X, U, P, h)
%OC_BACKWARD RK4 adjoint sweep backwards from zeta(T) = 0.
n = size(U,1) - 1; Z = zeros(n+1, 6);
for k = n+1:-1:2
    xm = 0.5*(X(k,:) + X(k-1,:))'; um = 0.5*(U(k,:) + U(k-1,:))'; z = Z(k,:)';
    k1 = oc_adj(X(k,:)', U(k,:)', z, P); k2 = oc_adj(xm, um, z - h/2*k1, P);
    k3 = oc_adj(xm, um, z - h/2*k2, P); k4 = oc_adj(X(k-1,:)', U(k-1,:)', z - h*k3, P);
    Z(k-1,:) = (z - h/6*(k1 + 2*k2 + 2*k3 + k4))';
end
end
