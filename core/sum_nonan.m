function s = sum_nonan(v)
%SUM_NONAN Sum ignoring NaN (works in MATLAB and Octave).
s = sum(v(~isnan(v)));
end
