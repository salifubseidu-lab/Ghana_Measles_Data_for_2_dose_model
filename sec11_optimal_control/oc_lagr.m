function L = oc_lagr(X, U, P)
%OC_LAGR Running cost of the objective at each grid point (rows of X, U).
I = X(:,5); extra = U(:,3).*P.alpha*P.K.*I.^2./(1 + P.K*I);
L = P.w*P.eps*X(:,4) + P.c1*U(:,1).*X(:,1) + P.c2*U(:,2).*X(:,2) + P.c3*extra ...
    + 0.5*(P.b1*U(:,1).^2 + P.b2*U(:,2).^2 + P.b3*U(:,3).^2);
end
