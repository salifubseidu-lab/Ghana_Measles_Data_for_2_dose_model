%S05_VERIFY_CLASSIFICATION  Random-parameter verification of Theorems 3.3, 5.1, 5.3,
%   Corollary 6.7 and the backward-bifurcation criterion K > Kc (Appendix A).
%   Set NS (number of samples) and MODE ('measles' or 'wide') before running.
setup_paths();
if ~exist('NS', 'var'), NS = 2000; end
if ~exist('MODE', 'var'), MODE = 'measles'; end
rng(1);
st = struct('n',0,'ngm_err',0,'dfe_ok',0,'Gam_min',Inf,'three',0,'two_Rcgt1',0, ...
    'two_KleKc',0,'upper_unstable',0,'lower_stable',0,'hopf_pairs',0,'newton_mismatch',0, ...
    'bb_tested',0,'bb_ok',0,'h2_fail',0,'detsign_tested',0,'detsign_ok',0,'conc_viol',0);
for it = 1:NS
    P = random_params(MODE); st.n = st.n + 1; rc = Rc_fun(P);
    st.ngm_err = max(st.ngm_err, abs(ngm_radius(P) - rc)/rc);
    st.dfe_ok = st.dfe_ok + ((max(real(eig(model_jac(dfe_state(P), P)))) < 0) == (rc < 1));
    st.Gam_min = min(st.Gam_min, Gamma_fun(P));
    lams = endemic_lambdas(P); n = numel(lams);
    if n >= 3, st.three = st.three + 1; end
    if n == 2 && rc > 1, st.two_Rcgt1 = st.two_Rcgt1 + 1; end
    if n == 2 && P.K <= Kc_fun(P), st.two_KleKc = st.two_KleKc + 1; end
    for j = 1:n
        X = ee_from_lambda(lams(j), P);
        assert(min(X) > 0 && max(abs(model_rhs(X, P))) < 1e-6*P.Lambda);
        ev = eig(model_jac(X, P));
        if j == n && max(real(ev)) > 0, st.upper_unstable = st.upper_unstable + 1;
            if sum(real(ev) > 0) >= 2 && any(abs(imag(ev(real(ev) > 0))) > 1e-9), st.hopf_pairs = st.hopf_pairs + 1; end
        end
        if j < n && max(real(ev)) < 0, st.lower_stable = st.lower_stable + 1; end
        % sign(det J) = -sign h'(lambda*)  (Theorem 6.6 / Corollary 6.7)
        hp = (h_fun(lams(j)*(1+1e-6), P) - h_fun(lams(j)*(1-1e-6), P))/(2e-6*lams(j));
        if abs(hp) > 1e-8
            st.detsign_tested = st.detsign_tested + 1;
            st.detsign_ok = st.detsign_ok + (sign(det(model_jac(X, P))) == -sign(hp));
        end
    end
    % concavity of h (Lemma 5.2(4)): analytic h'' (Theorem 5.3, Step 2) checked against
    % a finite-difference second derivative, and its sign
    lg = P.beta*P.mu*logspace(-4, 1, 7)/max(P.eps,1);
    for l = lg
        dl = 1e-3*l; hv = h_fun([l-dl, l, l+dl], P); hn = (hv(3) - 2*hv(2) + hv(1))/dl^2; ha = h2_fun(l, P);
        tolp = 1e-3*abs(ha) + 1e3*eps*max(abs(hv))/dl^2;      % truncation + round-off floor
        st.h2_fail = st.h2_fail + (abs(hn - ha) > tolp);
        if ha > 0, st.conc_viol = st.conc_viol + 1; end
    end
    % independent Newton check of equilibrium count (first 400 sets)
    if it <= 400
        x0 = dfe_state(P); N0 = P.Lambda/P.mu; found = [];
        for tr = 1:30
            Ig = N0*10^(-7 + 6.5*rand);
            g = [x0(1)*(0.01 + 0.99*rand); x0(2)*(0.01 + 0.99*rand); x0(3)*(0.1 + 0.9*rand); 0.8*Ig; Ig; N0*0.9*rand];
            [z, ok] = newton_ee(max(g, 1e-12), P, 80);
            if ok && min(z) > -1e-8*N0 && z(5) > 1e-9*N0
                if isempty(found) || all(abs(z(5) - found)./found > 1e-5), found(end+1) = z(5); end %#ok<SAGROW>
            end
        end
        Ic = arrayfun(@(l) ee_from_lambda(l, P)'*[0;0;0;0;1;0], lams);
        for f = found, if isempty(Ic) || all(abs(f - Ic)./Ic > 1e-4), st.newton_mismatch = st.newton_mismatch + 1; break; end, end
    end
    % backward bifurcation criterion: EE exist just below Rc = 1 iff K > Kc
    if P.K > 0 && abs(P.K/Kc_fun(P) - 1) > 0.1
        Q = P; Q.beta = beta_for_Rc(P, 1 - 1e-7); has = ~isempty(endemic_lambdas(Q));
        st.bb_tested = st.bb_tested + 1; st.bb_ok = st.bb_ok + (has == (P.K > Kc_fun(P)));
    end
end
disp(st);
save(fullfile(results_dir(), sprintf('verify_classification_%s.mat', MODE)), 'st');
