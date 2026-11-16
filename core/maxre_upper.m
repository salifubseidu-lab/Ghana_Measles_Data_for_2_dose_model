function m = maxre_upper(P)
%MAXRE_UPPER Largest real part of the Jacobian eigenvalues at the upper endemic equilibrium.
[~, ev] = upper_ee(P);
if isempty(ev), m = NaN; else, m = max(real(ev)); end
end
