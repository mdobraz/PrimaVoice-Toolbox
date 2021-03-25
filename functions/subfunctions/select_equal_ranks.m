function selection = select_equal_ranks(values_or_max_rank,n_groups,iterations)
if ~exist('iterations','var')
    iterations = 1000;
end

if length(values_or_max_rank) > 1
    values = values_or_max_rank;
    max_rank = length(values);
elseif length(values_or_max_rank) == 1
    max_rank = values_or_max_rank;
else
    error('''values_or_max_rank'' must be of length=1 to specify a max rank or length>1 to specify values')
end

n_elements = floor(max_rank / n_groups);

all_sel = nan(n_elements,n_groups,iterations);
diffs = nan(iterations,1);
for i = 1:iterations
    all_sel(:,:,i) = random_selection(max_rank,n_groups,n_elements);
    if exist('values','var')
        diffs(i) = sum(abs(diff(sum(values(all_sel(:,:,i)))))); % minimize value difference between groups
    else
        diffs(i) = sum(abs(diff(sum(all_sel(:,:,i))))); % minimize rank difference between groups
    end
end

[~,I] = min(diffs);
selection = sort(all_sel(:,:,I));




function rand_sel = random_selection(max_rank,n_groups,n_elements)
rand_sel = nan(n_elements,n_groups);
remaining = (1:max_rank)';
candidat = nan;
for e = 1:n_elements
    for g = 1:n_groups        
        remaining(remaining == candidat) = [];
        candidat = remaining(randi(length(remaining)));
        rand_sel(e,g) = candidat;
    end
end

