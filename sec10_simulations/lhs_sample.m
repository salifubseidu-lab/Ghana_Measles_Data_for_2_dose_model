function U = lhs_sample(n, d)
%LHS_SAMPLE Latin hypercube sample of n points in [0,1]^d (no toolbox needed).
U = zeros(n, d);
for j = 1:d, U(:,j) = (randperm(n)' - rand(n,1))/n; end
end
