function [xbest, fbest, info] = de_optimize(fun, bounds, opts)
%DE_OPTIMIZE Differential evolution (rand/best/1/bin, dithered mutation), minimising
%   fun, which maps a d-by-m matrix of candidates to a 1-by-m vector (vectorised).
%   opts: popsize (25, times d), maxiter (600), F ([0.5 1]), CR (0.7), tol (1e-8), seed.
d = size(bounds,1); lo = bounds(:,1); hi = bounds(:,2);
if ~isfield(opts,'popsize'), opts.popsize = 25; end
if ~isfield(opts,'maxiter'), opts.maxiter = 600; end
if ~isfield(opts,'F'), opts.F = [0.5 1]; end
if ~isfield(opts,'CR'), opts.CR = 0.7; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'verbose'), opts.verbose = false; end
if isfield(opts,'seed'), rng(opts.seed); end %#ok<RAND>
np = opts.popsize*d;
U = zeros(d, np); for j = 1:d, U(j,:) = (randperm(np) - rand(1,np))/np; end   % Latin hypercube start
X = lo*ones(1,np) + U.*((hi - lo)*ones(1,np)); f = fun(X);
for g = 1:opts.maxiter
    [fb, ib] = min(f); F = opts.F(1) + (opts.F(2) - opts.F(1))*rand;
    r1 = zeros(1,np); r2 = r1;
    for i = 1:np, rr = randperm(np, 3); rr(rr == i) = []; r1(i) = rr(1); r2(i) = rr(2); end
    V = X(:,ib)*ones(1,np) + F*(X(:,r1) - X(:,r2));                   % best/1
    cross = rand(d, np) < opts.CR; jr = randi(d, 1, np);
    cross(sub2ind([d np], jr, 1:np)) = true;
    T = X; T(cross) = V(cross);
    T = min(max(T, lo*ones(1,np)), hi*ones(1,np));
    ft = fun(T); better = ft <= f; X(:,better) = T(:,better); f(better) = ft(better);
    if opts.verbose && mod(g, 50) == 0, fprintf('  DE gen %d: best %.4f\n', g, min(f)); end
    if std(f) <= opts.tol*abs(mean(f)), break; end
end
[fbest, ib] = min(f); xbest = X(:,ib); info.gen = g; info.pop = X; info.f = f;
end
