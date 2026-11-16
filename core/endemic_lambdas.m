function lams = endemic_lambdas(P)
%ENDEMIC_LAMBDAS Positive real roots of the cubic (forces of infection of the
%   endemic equilibria), sorted ascending.
c = cubic_coeffs(P);
while ~isempty(c) && c(1) == 0, c = c(2:end); end
r = roots(c); lams = [];
for i = 1:numel(r)
    z = r(i);
    if abs(imag(z)) < 1e-9*max(1, abs(z)) && real(z) > 1e-12
        lams(end+1) = real(z); %#ok<AGROW>
    end
end
lams = sort(lams);
end
