function Us = oc_ustar(X, Z, P, umax)
%OC_USTAR Characterisation of the optimal controls (Theorem 11.4), projected on [0, umax].
I = X(:,5);
u1 = X(:,1).*(Z(:,1) - Z(:,2) - P.c1)/P.b1;
u2 = X(:,2).*(Z(:,2) - Z(:,3) - P.c2)/P.b2;
u3 = P.alpha*P.K*I.^2.*(Z(:,5) - Z(:,6) - P.c3)./(P.b3*(1 + P.K*I));
Us = min(max([u1 u2 u3], 0), ones(size(X,1),1)*umax(:)');
end
