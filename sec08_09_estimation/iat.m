function tau = iat(x)
%IAT Integrated autocorrelation time of a chain averaged over walkers (Sokal window,
%   c = 5, as in emcee).  x: nsteps-by-nw.
[n, nw] = size(x); f = zeros(n, 1);
for w = 1:nw
    y = x(:,w) - mean(x(:,w)); nf = 2^nextpow2(2*n);
    F = fft(y, nf); ac = real(ifft(F.*conj(F))); ac = ac(1:n)/ac(1); f = f + ac;
end
f = f/nw; taus = 2*cumsum(f) - 1; win = find((1:n)' >= 5*taus, 1);
if isempty(win), win = n; end
tau = taus(win);
end
