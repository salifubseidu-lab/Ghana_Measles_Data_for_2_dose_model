%S06_HOPF_TABLE  Table 8: Hopf thresholds K_H1, K_H2 for several Rc (beta varied),
%   with 1/K_H and the period 2*pi/omega of the cycle born there; first Lyapunov
%   coefficients and the normal-form amplitude at Rc = 3.92 (Finding 6.10).
setup_paths();
P0 = baseline_params(); kc = Kc_fun(P0);
rcs = [1.2 1.5 2 3 3.92 5 8]; T = zeros(numel(rcs), 7);
for r = 1:numel(rcs)
    P = P0; P.beta = beta_for_Rc(P, rcs(r));
    Ks = kc*logspace(-4, 2, 121); v = zeros(size(Ks));
    for i = 1:numel(Ks), Q = P; Q.K = Ks(i); v(i) = maxre_upper(Q); end
    idx = find(diff(sign(v)) ~= 0); row = [rcs(r) NaN NaN NaN NaN NaN NaN];
    for j = 1:min(2, numel(idx))
        f = @(K) maxre_upper(setfield(P, 'K', K)); %#ok<SFLD>
        KH = fzero(f, [Ks(idx(j)) Ks(idx(j)+1)], optimset('TolX', 1e-14));
        Q = P; Q.K = KH; [~, ev] = upper_ee(Q); [~, ii] = max(real(ev)); om = abs(imag(ev(ii)));
        row(2 + 3*(j-1):4 + 3*(j-1)) = [KH/kc, 1/KH, 2*pi/om];
    end
    T(r,:) = row;
    fprintf('Rc=%5.2f  KH1/Kc=%.4f  1/KH1=%6.0f  period=%5.2f | KH2/Kc=%6.2f  1/KH2=%6.2f  period=%5.2f\n', T(r,:));
end
dlmwrite(fullfile(results_dir(), 'tab08_hopf_thresholds.csv'), T, 'precision', '%.6g');
% Lyapunov coefficients and normal-form amplitude at Rc = 3.92
P = P0; P.beta = beta_for_Rc(P, 3.92); r = find(rcs == 3.92);
for j = 1:2
    KH = T(r, 2 + 3*(j-1))*kc; Q = P; Q.K = KH; X = upper_ee(Q); [L1, om, ~, q] = lyap1(X, Q);
    dK = KH*1e-5; s = (maxre_upper(setfield(Q,'K',KH+dK)) - maxre_upper(setfield(Q,'K',KH-dK)))/(2*dK); %#ok<SFLD>
    fprintf('K_H%d: l1 = %.4g, omega = %.4g, dRe/dK = %.3g\n', j, L1, om, s);
    if j == 1
        sig = s*KH*0.01; amp = 2*abs(q(5))*sqrt(max(-sig/(om*L1), 0));
        fprintf('  normal-form I-amplitude at K = 1.01 K_H1: %.4g\n', amp);
        KH1 = KH; l1_KH1 = L1; amp_nf = amp; %#ok<NASGU>
    else
        l1_KH2 = L1; %#ok<NASGU>
    end
end
save(fullfile(results_dir(), 'hopf.mat'), 'T', 'KH1', 'l1_KH1', 'l1_KH2', 'amp_nf');
