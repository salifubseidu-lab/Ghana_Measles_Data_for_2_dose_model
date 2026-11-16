function dz = fit_rhs_cont(t, z, D, p, nu1, nu2, beta0, M, k40)
%FIT_RHS_CONT Right-hand side of the fitting model with cumulative incidence (for the
%   solver check), drivers piecewise constant per day.
i = min(floor(t*365.25 + 1e-9) + 1, D.n_days - 1);
b = beta0*(1 + p.a*cos(2*pi*(t - p.phi))); if t >= 8.25 && t < 10, b = b*(1 - p.c); end
S = z(1); V1 = z(2); V2 = z(3); E = z(4); I = z(5); R = z(6); N = S + V1 + V2 + E + I + R;
lam = b*(I + p.iota)/N; mu = D.mu(i); th = p.theta; ell = M.ell;
dz = [D.Lam(i) + ell*V1 - lam*S - (nu1(i) + mu)*S;
      nu1(i)*S - th*lam*V1 - (ell + nu2(i) + mu)*V1;
      nu2(i)*V1 - mu*V2;
      lam*(S + th*V1) - (M.EPS + mu)*E;
      M.EPS*E - (k40 + mu)*I;
      (k40 - M.DD)*I - mu*R;
      M.EPS*E];
end
