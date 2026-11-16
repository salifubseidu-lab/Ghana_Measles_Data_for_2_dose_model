function [r, p] = prcc_calc(X, y)
%PRCC_CALC Partial rank correlation coefficients of output y with each column of X,
%   with two-sided t-test p-values (df = n - 2 - (k-1)).  No toolbox needed.
n = size(X,1); k = size(X,2);
R = zeros(n,k); for j = 1:k, R(:,j) = rank_vec(X(:,j)); end
ry = rank_vec(y); r = zeros(1,k); p = r;
for j = 1:k
    Z = [ones(n,1), R(:, [1:j-1, j+1:k])];
    ex = R(:,j) - Z*(Z\R(:,j)); ey = ry - Z*(Z\ry);
    c = corrcoef(ex, ey); r(j) = c(1,2);
    df = n - 2 - (k - 1); t = r(j)*sqrt(df/(1 - r(j)^2));
    p(j) = betainc(df/(df + t^2), df/2, 0.5);            % = 2*(1 - tcdf(|t|, df))
end
end
function rk = rank_vec(v)
[~, idx] = sort(v); rk = zeros(size(v)); rk(idx) = (1:numel(v))';
end
