function y = nb_rand(k, m)
%NB_RAND Negative-binomial draws with size k and mean m (elementwise), generated as a
%   gamma-Poisson mixture.  Core MATLAB only (no Statistics Toolbox: own gamma and
%   Poisson samplers).
lam = gamma_rand(k.*ones(size(m))).*(m./k);
y = pois_rand(lam);
end
function g = gamma_rand(a)
% Marsaglia-Tsang (2000) sampler for Gamma(a, 1); a < 1 via the boost a+1 and U^(1/a)
g = zeros(size(a));
for i = 1:numel(a)
    ai = a(i); boost = 1;
    if ai < 1, boost = rand^(1/ai); ai = ai + 1; end
    d = ai - 1/3; c = 1/sqrt(9*d);
    while true
        x = randn; v = (1 + c*x)^3;
        if v <= 0, continue; end
        u = rand;
        if log(u) < 0.5*x^2 + d - d*v + d*log(v), break; end
    end
    g(i) = d*v*boost;
end
end
function y = pois_rand(lam)
y = zeros(size(lam));
small = lam < 30;
for i = find(small(:))'                         % multiplication method
    L = exp(-lam(i)); kk = 0; p = rand;
    while p > L, kk = kk + 1; p = p*rand; end
    y(i) = kk;
end
for i = find(~small(:))'                        % PTRS (Hormann 1993)
    l = lam(i); slam = sqrt(l); loglam = log(l); b = 0.931 + 2.53*slam; a = -0.059 + 0.02483*b;
    invalpha = 1.1239 + 1.1328/(b - 3.4); vr = 0.9277 - 3.6224/(b - 2);
    while true
        U = rand - 0.5; V = rand; us = 0.5 - abs(U); kk = floor((2*a/us + b)*U + l + 0.43);
        if us >= 0.07 && V <= vr, break; end
        if kk < 0 || (us < 0.013 && V > us), continue; end
        if log(V) + log(invalpha) - log(a/(us*us) + b) <= -l + kk*loglam - gammaln(kk + 1), break; end
    end
    y(i) = kk;
end
end
