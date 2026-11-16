function p = est_unpack(M, TH)
%EST_UNPACK Natural-scale parameters from the columns of TH (npar-by-m).
p.R0 = TH(1,:); p.a = TH(2,:); p.phi = TH(3,:); p.iota = 10.^TH(4,:); p.rho = 1./(1 + exp(-TH(5,:)));
p.s0 = TH(6,:); p.w1 = TH(7,:); p.k = exp(TH(8,:)); p.theta = 1 - p.w1;
if isempty(M.kappa_fixed), p.kappa = TH(9,:); else, p.kappa = M.kappa_fixed*ones(1, size(TH,2)); end
if M.covid, p.c = TH(10,:); else, p.c = zeros(1, size(TH,2)); end
if M.sia_est, p.h = TH(end,:); else, p.h = ones(1, size(TH,2)); end
end
