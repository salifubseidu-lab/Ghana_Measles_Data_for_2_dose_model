function Q = P_on_branch(lam, P)
%P_ON_BRANCH Parameter set whose equilibrium branch passes through lam.
Q = P; Q.beta = qhat_branch(lam, P)/P.mu;
end
