function dy = ode_log(~, y, P)
%ODE_LOG Model in logarithmic coordinates y = log(x) (for stiff, deep troughs).
x = exp(max(y, -700));          % clamp: avoids exp underflow to 0 (I < 1e-304)
dy = model_rhs(x, P)./x;
end
