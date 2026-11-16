function v = prctile_(x, pct)
%PRCTILE_ Percentiles along columns (linear interpolation; no toolbox needed).
%   For a row vector x the percentiles of x are returned.
if isrow(x), x = x(:); end
x = sort(x, 1); n = size(x, 1); v = zeros(numel(pct), size(x, 2));
for i = 1:numel(pct)
    r = 1 + (n - 1)*pct(i)/100; lo = floor(r); hi = ceil(r); w = r - lo;
    v(i,:) = (1 - w)*x(lo,:) + w*x(hi,:);
end
if size(x,2) == 1, v = v.'; end
end
