function L = first_l1(P0, K)
%FIRST_L1 l1 at the first Hopf point on the branch for fixed K (used to locate GH points).
P = P0; P.K = K; H = hopf_points(P, 1.5, 800);
if isempty(H), L = NaN; else, L = H(1,4); end
end
