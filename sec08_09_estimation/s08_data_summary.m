%S08_DATA_SUMMARY  Section 8: annual surveillance table (Table 10), demographic and
%   coverage drivers, age and vaccination status of 2022-24 cases, mean age (prior for
%   lambda_ref) and the screening-method vaccine effectiveness.
setup_paths();
A = read_csv_table(fullfile(data_dir(), 'ghana_annual_surveillance.csv'));
fprintf('Year  Susp   Lab  Epi  Clin  Total  Inc/1e6  Discarded/1e5\n');
for i = 1:numel(A.year)
    fprintf('%d %5d %5d %4d %5d %6d %8.2f %8.2f\n', A.year(i), A.suspected(i), A.lab(i), A.epi(i), A.clin(i), A.total(i), A.incidence_per_million(i), A.discarded_rate_per_100k(i));
end
D = est_inputs('all'); a = D.annual;
fprintf('\nYear  population   CBR(per mil)  growth   mu      p1     p2c\n');
for i = 1:numel(a.year), fprintf('%d %12.0f %8.2f %8.4f %7.4f %6.3f %6.3f\n', a.year(i), a.pop(i), 1000*a.cbr(i), a.growth(i), a.mu(i), a.p1(i), a.p2c(i)); end
fprintf('fit months %d, hold-out months %d (missing months excluded)\n', sum(D.fit & ~isnan(D.y)), sum(~D.fit & ~isnan(D.y)));
% outbreak table: age and vaccination status, 2022-24
O = read_csv_table(fullfile(data_dir(), 'ghana_outbreak_age_vaccination.csv'));
fn = fieldnames(O); sel = strcmp(O.CODE, 'GHA') & O.YEAR >= 2022 & O.YEAR <= 2024;
age = {'lt_1_year','x1_to_4_years','x5_to_9_years','x10_to_14_years','gteq_15_years'}; mid = [0.5 2.5 7 12 22];
have = all(isfield(O, age));
if have
    n = cellfun(@(c) sum_nonan(O.(c)(sel)), age);
    fprintf('\ncases by age 2022-24: %s (n = %d), shares %s\n', mat2str(n), sum(n), mat2str(round(1000*n/sum(n))/10));
    fprintf('mean age (interval midpoints, 22 y for 15+): %.2f y -> lambda_ref prior median 1/%.1f\n', sum(n.*mid)/sum(n), sum(n.*mid)/sum(n));
else
    fprintf('\ncolumns of the outbreak table: %s\n', strjoin(fn', ', '));
end
vac = {'vaccinated','unvaccinated'};
if all(isfield(O, vac))
    v = sum_nonan(O.vaccinated(sel)); u = sum_nonan(O.unvaccinated(sel)); un = sum_nonan(O.unknown(sel));
    fprintf('vaccinated %d, unvaccinated %d, unknown %d\n', v, u, un); PCV = v/(v + u);
    for PPV = [0.92 0.87]
        fprintf('PCV = %.3f, PPV = %.2f: screening VE = %.3f\n', PCV, PPV, 1 - PCV/(1 - PCV)*(1 - PPV)/PPV);
    end
end
