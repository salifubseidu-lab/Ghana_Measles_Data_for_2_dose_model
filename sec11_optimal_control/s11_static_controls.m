%S11_STATIC_CONTROLS  Proposition 11.1 and Finding 11.2 (baseline, K = 2 Kc, Rc = 3.92):
%   constant u3 maps (gamma, alpha) -> (gamma + alpha u3, alpha (1-u3)); Kc(u3) = Kc/(1-u3).
%   Thresholds: u3 removing bistability (R_SN -> 1) and u3 restabilising the endemic state.
setup_paths();
P = oc_params(); kc = Kc_fun(P);
mk = @(u3) setfield(setfield(P, 'gamma', P.gamma + P.alpha*u3), 'alpha', P.alpha*(1 - u3)); %#ok<SFLD>
for u3 = [0 0.25 0.5]
    Q = mk(u3); fprintf('u3 = %.2f: Kc(u3) = %.4e (Kc/(1-u3) = %.4e), Rc = %.4f, R_SN = %.4f\n', ...
        u3, Kc_fun(Q), kc/(1 - u3), Rc_fun(Q), RSN_fun(Q));
end
fprintf('bistability removed for u3 >= 1 - Kc/K = %.3f\n', 1 - kc/P.K);
g = @(u3) maxre_upper(mk(u3));
u3s = fzero(g, [0 0.95]);
fprintf('endemic equilibrium restabilised for u3 >= %.3f (Rc stays %.4f)\n', u3s, Rc_fun(mk(u3s)));
