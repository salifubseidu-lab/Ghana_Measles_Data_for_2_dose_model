function M = est_setup(varargin)
%EST_SETUP Estimation model of Section 9 (options as name/value pairs):
%   'target' ('all'|'conf'), 'cov' ('wuenic'|'dhs'), 'ell' (0.01), 'eta_scale' (1),
%   'covid' (true), 'sia_est' (true), 'kappa_fixed' ([] = estimated),
%   'prior_R0' ([15 3]), 'prior_omega' ([0.80 0.06]), 'obj' ('post'|'lik').
o = struct('target','all','cov','wuenic','ell',0.01,'eta_scale',1,'covid',true,'sia_est',true, ...
           'kappa_fixed',[],'prior_R0',[15 3],'prior_omega',[0.80 0.06],'obj','post');
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
M = o; M.D = est_inputs(o.target, o.cov, o.eta_scale);
M.YD = 365.25; M.EPS = M.YD/10; M.K40T = M.YD/7; M.CFR = 0.02; M.DD = M.CFR/(1 - M.CFR)*M.K40T;
M.KAPPA_REF = 1/5.2; M.RHO_PRIOR = [log(0.14/0.86) 1.0];
y = M.D.y; M.mask = M.D.fit & ~isnan(y); M.y = y; M.y(isnan(y)) = 0;
lo = [6 0 0 -2 -9 0.002 0.5 log(0.2) 0.0]; hi = [30 0.6 1 3 4 0.30 0.99 log(200) 1.5];
names = {'R0','a','phi','log10iota','logitrho','sigma0','omega1','logk','kappa'};
if o.covid, lo(end+1) = 0; hi(end+1) = 0.9; names{end+1} = 'c'; end
if o.sia_est, lo(end+1) = 0.05; hi(end+1) = 1.6; names{end+1} = 'h'; end
M.bounds = [lo' hi']; M.names = names; M.npar = numel(lo);
end
