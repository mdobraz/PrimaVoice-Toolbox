function X = clean_binary(X,min_n_consec_ones,min_n_consec_zeros)
% removes groups of consecutive ones or zeros in X, whose size is below
% 'min_n_consec_ones' or 'min_n_consec_zeros'

% the algorithm is conservative, it will not change small groups of zeros
% to ones if these groups follow or are followed by small a small group of
% one.

% calculate group sizes
group_sizes = get_group_sizes(X);

to_remove = (group_sizes(:,1) < min_n_consec_ones & group_sizes(:,2) == 1) | (group_sizes(:,1) < min_n_consec_zeros & group_sizes(:,2) == 0);

for i = 1 :length(to_remove)
    if to_remove(i)
        if group_sizes(i,2) % if it's a group of ones
            idx1 = sum(group_sizes(1:i)) - group_sizes(i) + 1;
            idx2 = idx1 + group_sizes(i) - 1;
            X(idx1:idx2) = 0; % replace group by zeros
        else % if it's a group of zeros
            ok = 0;
            if i > 1
                if ~to_remove(i-1) % previous group was not removed
                    ok = 1;
                end
            end
            if i < length(to_remove)
                if ~to_remove(i+1) && ok % next group will not be removed and previous was not removed
                    ok = 1;
                else
                    ok = 0;
                end
            end
            if ok % if previous and next group of ones has not to be removed
                idx1 = sum(group_sizes(1:i)) - group_sizes(i) + 1;
                idx2 = idx1 + group_sizes(i) - 1;
                X(idx1:idx2) = 1; % replace group by ones
            end
        end
    end
end

% i = 1;
% while i <= length(X)
%     if X(i)
%         n_ones = find(~X(i:end),1) - 1;
%         if isempty(n_ones)
%             n_ones = length(X) - i + 1;
%         end
%         if n_ones < min_n_consec_ones
%             end_idx = min([i+n_ones length(X)]);
%             X(i:end_idx) = 0; % remove small groups of ones 
%         end
%         i = i + n_ones;
%     else
%         i = i + 1;
%     end
% end