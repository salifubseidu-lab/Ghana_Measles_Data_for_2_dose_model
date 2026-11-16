function [inc, st] = est_simulate(M, TH, sia_frac)
%EST_SIMULATE Fitting model (9.1) for every column of TH (vectorised RK4, step 1 day).
%   inc: monthly incidence of infectiousness C_m (n_months-by-m);
%   st:  states at month starts (n_starts-by-6-by-m).
D = M.D; p = est_unpack(M, TH); m = size(TH,2);
if nargin < 3 || isempty(sia_frac), sia_frac = D.sia_frac; end
YD = M.YD; dt = 1/YD; EPS = M.EPS; DD = M.DD; k40 = M.K40T + DD; ell = M.ell;
% vaccination rates by the competing-hazard mapping (Lemma 8.2), per day and parameter set
nu1 = (D.p1.*(D.mu))*ones(1,m) + D.p1*p.kappa;   nu1 = nu1./((1 - D.p1)*ones(1,m));
nu2 = (D.p2.*(ell + D.mu))*ones(1,m) + D.p2*(p.theta.*p.kappa); nu2 = nu2./((1 - D.p2)*ones(1,m));
% initial state (Section 9.2)
mu0 = D.mu(1); n1 = nu1(1,:); n2 = nu2(1,:); k2 = ell + n2 + mu0; Dd = mu0*(mu0 + ell + n1 + n2) + n1.*n2;
fS = mu0*k2./Dd; fV1 = mu0*n1./Dd; fV2 = n1.*n2./Dd; N0 = D.N(1);
S = p.s0*N0./(1 + p.theta.*fV1./fS); V1 = S.*fV1./fS; V2 = min(fV2*N0, N0 - S - V1);
y1 = max(M.y(1), 1); E = 12*y1./(p.rho*EPS); I = E*EPS/(M.K40T + DD); R = N0 - S - V1 - V2 - E - I;
x = [S; V1; V2; E; I; R];
beta0 = p.R0*(EPS + mu0)*(k40 + mu0)/EPS;
sf = min(sia_frac(:)*p.h, 0.99);                                    % pulse fractions (n_sia-by-m)
n = D.n_days; daily = zeros(n, m); ms = D.mstart; st = zeros(numel(ms), 6, m);
mi = 1; si = 1; covlo = 2020.25 - 2012; covhi = 2022 - 2012;
for i = 1:n-1
    if mi <= numel(ms) && i == ms(mi), st(mi,:,:) = reshape(x, 1, 6, m); mi = mi + 1; end
    if si <= numel(D.sia_idx) && i == D.sia_idx(si)
        dS = sf(si,:).*x(1,:); dV = sf(si,:).*x(2,:);
        x(1,:) = x(1,:) - dS; x(2,:) = x(2,:) + dS - dV; x(3,:) = x(3,:) + dV; si = si + 1;
    end
    tt = D.t(i);
    b = beta0.*(1 + p.a.*cos(2*pi*(tt - p.phi)));
    if tt >= covlo && tt < covhi, b = b.*(1 - p.c); end
    L = D.Lam(i); mu = D.mu(i); v1 = nu1(i,:); v2 = nu2(i,:);
    [k1, f1] = rhs_fit(x, b, p.theta, L, mu, v1, v2, ell, k40, p.iota, EPS, DD);
    [k2v, f2] = rhs_fit(x + 0.5*dt*k1, b, p.theta, L, mu, v1, v2, ell, k40, p.iota, EPS, DD);
    [k3, f3] = rhs_fit(x + 0.5*dt*k2v, b, p.theta, L, mu, v1, v2, ell, k40, p.iota, EPS, DD);
    [k4, f4] = rhs_fit(x + dt*k3, b, p.theta, L, mu, v1, v2, ell, k40, p.iota, EPS, DD);
    x = max(x + dt/6*(k1 + 2*k2v + 2*k3 + k4), 0);
    daily(i,:) = dt/6*(f1 + 2*f2 + 2*f3 + f4);
end
nm = numel(ms) - 1; inc = zeros(nm, m);
for k = 1:nm, inc(k,:) = sum(daily(ms(k):ms(k+1)-1, :), 1); end
end

function [dx, f] = rhs_fit(x, b, th, L, mu, v1, v2, ell, k40, iota, EPS, DD)
S = x(1,:); V1 = x(2,:); V2 = x(3,:); E = x(4,:); I = x(5,:); R = x(6,:);
N = S + V1 + V2 + E + I + R; lam = b.*(I + iota)./N;
dx = [L + ell*V1 - lam.*S - (v1 + mu).*S;
      v1.*S - th.*lam.*V1 - (ell + v2 + mu).*V1;
      v2.*V1 - mu*V2;
      lam.*(S + th.*V1) - (EPS + mu)*E;
      EPS*E - (k40 + mu)*I;
      (k40 - DD)*I - mu*R];
f = EPS*E;
end
