function [chain, lp, acc] = ensemble_mcmc(logp, X0, nsteps, pde)
%ENSEMBLE_MCMC Affine-invariant ensemble sampler (Goodman-Weare 2010; emcee scheme with
%   two complementary half-ensembles).  Moves: differential-evolution move (probability
%   pde, gamma = 2.38/sqrt(2d)) and the stretch move (a = 2) otherwise.
%   logp maps d-by-m to 1-by-m.  chain: nsteps-by-d-by-nw; lp: nsteps-by-nw.
if nargin < 4, pde = 0.8; end
[d, nw] = size(X0); X = X0; L = logp(X); chain = zeros(nsteps, d, nw); lp = zeros(nsteps, nw);
half = {1:floor(nw/2), floor(nw/2)+1:nw}; nacc = 0; g0 = 2.38/sqrt(2*d); a = 2;
for it = 1:nsteps
    for hh = 1:2
        S = half{hh}; C = half{3 - hh}; ns = numel(S); nc = numel(C);
        if rand < pde                                   % DE move (symmetric)
            Q = zeros(d, ns);
            for i = 1:ns, pr = C(randperm(nc, 2)); Q(:,i) = X(:,S(i)) + g0*(1 + 1e-5*randn)*(X(:,pr(1)) - X(:,pr(2))); end
            lf = zeros(1, ns);
        else                                            % stretch move
            z = ((a - 1)*rand(1, ns) + 1).^2/a; j = C(randi(nc, 1, ns));
            Q = X(:,j) + (ones(d,1)*z).*(X(:,S) - X(:,j)); lf = (d - 1)*log(z);
        end
        Lq = logp(Q); ok = log(rand(1, ns)) < lf + Lq - L(S); ok(~isfinite(Lq)) = false;
        X(:,S(ok)) = Q(:,ok); L(S(ok)) = Lq(ok); nacc = nacc + sum(ok);
    end
    chain(it,:,:) = reshape(X, 1, d, nw); lp(it,:) = L;
end
acc = nacc/(nsteps*nw);
end
