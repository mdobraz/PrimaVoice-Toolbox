function group_sizes = get_group_sizes(X)

i = 1;
g = 1;
group_sizes = nan(1,2);
while i <= length(X)
    n_group = find(X(i:end) ~= X(i),1) - 1;
    if isempty(n_group)
        n_group = length(X) - i + 1;
    end
    group_sizes(g,1) = n_group;
    group_sizes(g,2) = X(i);
    g = g + 1;
    i = i + n_group;
end
