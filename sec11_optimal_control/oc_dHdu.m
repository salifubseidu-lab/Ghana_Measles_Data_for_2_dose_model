function G = oc_dHdu(X, U, Z, P)
%OC_DHDU Gradient of the Hamiltonian with respect to the controls at each grid point.
I = X(:,5); ex = P.alpha*P.K*I.^2./(1 + P.K*I);
G = [P.c1*X(:,1) + P.b1*U(:,1) + (Z(:,2) - Z(:,1)).*X(:,1), ...
     P.c2*X(:,2) + P.b2*U(:,2) + (Z(:,3) - Z(:,2)).*X(:,2), ...
     P.c3*ex + P.b3*U(:,3) + ex.*(Z(:,6) - Z(:,5))];
end
