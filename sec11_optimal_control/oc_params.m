function P = oc_params(Kfac, rc)
%OC_PARAMS Baseline (Table 7) with K = Kfac*Kc and Rc = rc, plus the illustrative
%   economic parameters of Section 11.5 (US$): w, c1, c2, c3, b1, b2, b3.
if nargin < 1, Kfac = 2; end
if nargin < 2, rc = 3.92; end
P = baseline_params(); P.K = Kfac*Kc_fun(P); P.beta = beta_for_Rc(P, rc);
P.w = 150; P.c1 = 1.5; P.c2 = 1.5; P.c3 = 20; P.b1 = 2e4; P.b2 = 2e4; P.b3 = 1e6;
end
