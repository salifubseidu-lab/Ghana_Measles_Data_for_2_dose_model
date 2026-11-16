function ci = profile_ci(grid, d)
%PROFILE_CI Interpolated interval where the profile difference d is below 1.92
%   (NaN or the grid end if not crossed).
thr = 1.92; [dm, i0] = min(d); ci = [grid(1) grid(end)];
if dm > thr, ci = [NaN NaN]; return; end      % grid does not bracket the optimum
for i = i0:-1:2, if d(i-1) > thr, ci(1) = grid(i-1) + (thr - d(i-1))*(grid(i) - grid(i-1))/(d(i) - d(i-1)); break; end, end
for i = i0:numel(d)-1, if d(i+1) > thr, ci(2) = grid(i) + (thr - d(i))*(grid(i+1) - grid(i))/(d(i+1) - d(i)); break; end, end
end
