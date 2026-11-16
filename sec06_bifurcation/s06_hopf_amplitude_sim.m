%S06_HOPF_AMPLITUDE_SIM  Direct simulation of the cycle amplitude at K = 1.01 K_H1
%   (Rc = 3.92), to compare with the normal-form prediction (Finding 6.10).
%   Long run (TEND years); the amplitude is (max I - min I)/2 over the last 300 years.
setup_paths();
if ~exist('TEND', 'var'), TEND = 8000; end
S = load(fullfile(results_dir(), 'hopf.mat'));
P = baseline_params(); P.beta = beta_for_Rc(P, 3.92); P.K = 1.01*S.KH1;
X = upper_ee(P); X0 = X.*(1 + 0.002*[1;-1;0;1;1;0]);
tt = linspace(TEND - 300, TEND, 20000);
Yl = ode_long(P, log(X0), tt, 1e-10, 1e-10); I = exp(Yl(:,5));
fprintf('simulated I-amplitude at 1.01 K_H1: %.4g (normal form %.4g)\n', (max(I) - min(I))/2, S.amp_nf);
