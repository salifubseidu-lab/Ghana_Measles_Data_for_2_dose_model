function D = est_inputs(target, cov, eta_scale)
%EST_INPUTS Build the observation series and the daily drivers of Section 8:
%   monthly cases y_m (target 'all' = lab + epi + clinical, or 'conf' = lab + epi),
%   population P(t) (exponential interpolation of mid-year values), births
%   Lambda(t) = b(t) P(t), exit rate mu(t) = b(t) - g(t) (Lemma 8.1), coverage p1(t),
%   p2c(t) (WUENIC, or 'dhs' = WUENIC - 3.2 points), and the SIA pulse days.
if nargin < 1, target = 'all'; end
if nargin < 2, cov = 'wuenic'; end
if nargin < 3, eta_scale = 1; end
dd = data_dir();
% ---- monthly cases (2012-01 .. 2026-05)
C = read_csv_table(fullfile(dd, 'ghana_monthly_cases.csv'));
y = C.lab_confirmed + C.epi_linked; if strcmp(target, 'all'), y = y + C.clinical; end
D.year = C.year; D.month = C.month; D.y = y;
D.fit = (C.year < 2024) | (C.year == 2024 & C.month <= 10);       % fit window 2012-01 .. 2024-10
D.datenum = datenum(C.year, C.month, 1);
% ---- demography 2010..2027 with the extrapolation rules of the paper
W = read_csv_table(fullfile(dd, 'ghana_demography.csv')); yrs = (2010:2027)';
b = nan(size(yrs)); p = b;
for i = 1:numel(yrs)                     % exact year lookup (interp1 would propagate NaN neighbours)
    k = find(W.year == yrs(i));
    if ~isempty(k), b(i) = W.cbr_per1000(k)/1000; p(i) = W.population(k); end
end
sel = yrs >= 2019 & yrs <= 2024; cf = polyfit(yrs(sel), b(sel), 1);
b(isnan(b)) = polyval(cf, yrs(isnan(b)));
for i = 3:numel(yrs), if isnan(p(i)), p(i) = p(i-1)^2/p(i-2); end, end
g = [log(p(2:end)./p(1:end-1)); NaN]; g(end) = g(end-1);           % growth y -> y+1 (forward fill)
mu = b - g;
% ---- coverage 2010..2027 (forward fill), conditional second-dose uptake
V = read_csv_table(fullfile(dd, 'ghana_coverage_wuenic.csv'));
c1 = nan(size(yrs)); c2 = c1;
for i = 1:numel(yrs)
    k = find(V.year == yrs(i)); if ~isempty(k), c1(i) = V.MCV1(k); c2(i) = V.MCV2(k); end
end
for i = 2:numel(yrs), if isnan(c1(i)), c1(i) = c1(i-1); end, if isnan(c2(i)), c2(i) = c2(i-1); end, end
c1 = c1/100; c2 = c2/100;
if strcmp(cov, 'dhs'), c1 = c1 - 0.032; c2 = c2 - 0.032; end
p2c = min(c2./c1, 0.98);
D.annual = struct('year', yrs, 'pop', p, 'cbr', b, 'growth', g, 'mu', mu, 'p1', c1, 'p2c', p2c);
% ---- daily grid 2012-01-01 .. 2026-06-01
t0 = datenum(2012,1,1); days = (t0:datenum(2026,6,1))'; n = numel(days);
dv = datevec(days); dyear = dv(:,1);
D.t = (days - t0)/365.25;
D.N = exp(interp1(yrs + 0.5, log(p), 2012 + D.t, 'linear', 'extrap'));
iy = dyear - 2009;                                                   % index into yrs
D.mu = mu(iy); D.Lam = b(iy).*D.N; D.p1 = c1(iy); D.p2 = p2c(iy);
% month starts (1-based day indices), 2012-01 .. 2026-06
[Y, M] = meshgrid(2012:2026, 1:12); ms = datenum(Y(:), M(:), 1); ms = ms(ms <= datenum(2026,6,1));
D.mstart = ms - t0 + 1;
% SIA pulses: day index (1-based), nominal fraction eta*coverage (scaled by eta_scale)
S = read_csv_table(fullfile(dd, 'ghana_sia.csv'));
D.sia_idx = datenum(S.year, S.month, S.day) - t0 + 1;
D.sia_frac = min(S.survey_coverage.*S.eta*eta_scale, 0.99);
D.n_days = n; D.days = days;
end
