function c = crps_sample(s, y)
%CRPS_SAMPLE Continuous ranked probability score of observation y for predictive sample s.
s = s(:); c = mean(abs(s - y)) - 0.5*mean(mean(abs(s - s.')));
end
