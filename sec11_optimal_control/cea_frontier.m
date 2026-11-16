function T = cea_frontier(keys, eff, cost)
%CEA_FRONTIER Strong and extended dominance and ICERs on the efficient frontier.
%   keys: cell of strategy names (first = comparator with eff = cost = 0).
%   Returns a struct array with fields key, eff, cost, acer, status ('', 'D', 'ED'), icer.
n = numel(keys); [~, o] = sortrows([eff(:) cost(:)]); keys = keys(o); eff = eff(o); cost = cost(o);
status = repmat({''}, 1, n);
for i = 1:n
    for j = 1:n
        if j ~= i && eff(j) >= eff(i) && cost(j) <= cost(i) && (eff(j) > eff(i) || cost(j) < cost(i)), status{i} = 'D'; end
    end
end
front = find(strcmp(status, ''));
changed = true;
while changed
    changed = false;
    ic = diff(cost(front))./diff(eff(front));
    for i = 2:numel(ic)
        if ic(i) < ic(i-1), status{front(i)} = 'ED'; front(i) = []; changed = true; break; end
    end
end
icer = nan(1, n); for i = 2:numel(front), icer(front(i)) = (cost(front(i)) - cost(front(i-1)))/(eff(front(i)) - eff(front(i-1))); end
T = struct('key', keys, 'eff', num2cell(eff), 'cost', num2cell(cost), 'acer', num2cell(cost./eff), 'status', status, 'icer', num2cell(icer));
end
