function inc = niamey_sim(TH, W)
%NIAMEY_SIM Weekly expected reported cases of the closed SEIR model for each column of
%   TH = [R0; s0; logit rho; log I0; log k] (RK4, 28 substeps per week; vectorised).
%   Population 7.7e5; 1/eps = 10 d, 1/k40 = 7 d.
NPOP = 770000; EPSW = 7/10; K40W = 7/7; m = size(TH,2);
b = TH(1,:)*K40W; rho = 1./(1 + exp(-TH(3,:))); I0 = exp(TH(4,:));
S = TH(2,:)*NPOP; E = I0*K40W/EPSW; I = I0; dt = 1/28; inc = zeros(W, m);
for w = 1:W
    c = zeros(1, m);
    for s = 1:28
        [a1, c1] = f(S, E, I);
        [a2, c2] = f(S + dt/2*a1(1,:), E + dt/2*a1(2,:), I + dt/2*a1(3,:));
        [a3, c3] = f(S + dt/2*a2(1,:), E + dt/2*a2(2,:), I + dt/2*a2(3,:));
        [a4, c4] = f(S + dt*a3(1,:), E + dt*a3(2,:), I + dt*a3(3,:));
        S = S + dt/6*(a1(1,:) + 2*a2(1,:) + 2*a3(1,:) + a4(1,:));
        E = E + dt/6*(a1(2,:) + 2*a2(2,:) + 2*a3(2,:) + a4(2,:));
        I = I + dt/6*(a1(3,:) + 2*a2(3,:) + 2*a3(3,:) + a4(3,:));
        c = c + dt/6*(c1 + 2*c2 + 2*c3 + c4);
    end
    inc(w,:) = c;
end
inc = inc.*(ones(W,1)*rho);
    function [d, cc] = f(S, E, I)
        lam = b.*I/NPOP;
        d = [-lam.*S; lam.*S - EPSW*E; EPSW*E - K40W*I]; cc = EPSW*E;
    end
end
