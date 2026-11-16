function T = read_csv_table(fn)
%READ_CSV_TABLE Read a simple comma-separated file with a header row into a struct of
%   column vectors (numeric where possible; non-numeric columns are returned as cell arrays).
txt = fileread(fn); lines = regexp(txt, '\r?\n', 'split'); lines = lines(~cellfun(@isempty, lines));
hdr = strsplit(lines{1}, ',', 'CollapseDelimiters', false); nr = numel(lines) - 1; C = repmat({''}, nr, numel(hdr));
for i = 1:nr, parts = strsplit(lines{i+1}, ',', 'CollapseDelimiters', false); C(i, 1:numel(parts)) = parts; end
for j = 1:numel(hdr)
    name = matlab.lang.makeValidName(strtrim(hdr{j}));
    v = str2double(C(:,j));
    if all(isnan(v) & ~strcmpi(strtrim(C(:,j)), 'nan') & ~cellfun(@isempty, C(:,j))), T.(name) = C(:,j);
    else, T.(name) = v; end
end
end
