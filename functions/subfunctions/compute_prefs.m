function [Prefs,Probs] = compute_prefs(conds,contrasts,dims,paths)

%% Get betas from vs_silence contrasts
nconds = length(conds);
all_ef = nan([dims nconds]);
all_ef_vect_cell = cell(nconds,1);
for i = 1:nconds
    % load ef_file
    ef_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{conds(i)} '_ef.nii.gz']);
    P = spm_vol(ef_file);
    Y = spm_read_vols(P); % put all contrasts in one variable
    % if MION, inverse values
    if strcmpi(paths.contrast_agent,'MION')
        Y = -Y;
    end
    all_ef(:,:,:,i) = Y;
    all_ef_vect_cell{i} = Y(logical(Y));
end

n_values = length(all_ef_vect_cell{1});
all_ef_vect = zeros(nconds,n_values);
for i = 1:nconds
    all_ef_vect(i,:) = all_ef_vect_cell{i};
end



%% Compute bootstrap distribution
max_bytes = 100000000; % max size of the difference vector: 100MB
n_perms = floor(max_bytes / n_values / 4);
n_diffs = n_perms * n_values;
fprintf('Computing bootstrap distribution based on %i permutations of each set of voxels, \nfor a total of %i differences...',n_perms,n_diffs)
D = single(nan(n_diffs,1));
i1 = 1;
i2 = n_values;
for k = 1:n_perms
    p = nan(nconds,n_values);
    for c = 1:nconds
        p(c,:) = randperm(n_values);
    end
    aev = nan(size(all_ef_vect));
    for c = 1:nconds
        aev(c,:) = all_ef_vect(c,p(c,:));
    end
    [a,ind] = sort(aev,'descend');

    D(i1:i2) = aev(1,:) - aev(2,:); % differences
    i1 = i2 + 1;
    i2 = i2 + n_values;
end






% %% Compute bootstrap distribution (old way)
% max_bytes = 100000000; % max size of the difference vector: 100MB
% cond_triang = ~diag(conds>0); % all possible differences between conditions
% n_cond_diffs = sum(cond_triang(:));
% n_perms = floor(max_bytes / n_values / n_cond_diffs / 4);
% n_diffs = n_perms * n_values * n_cond_diffs;
% fprintf('Computing bootstrap distribution based on %i permutations of each pair of conditions, \nfor a total of %i differences...',n_perms,n_diffs)
% D = single(nan(n_diffs,1));
% i1 = 1;
% i2 = n_values;
% for c1 = 1:nconds
%     for c2 = 1:nconds
%         if cond_triang(c1,c2)
%             for k = 1:n_perms
%                 a = 0;
%                 while any(a == 0) % exclude differences of the same indexes
%                     p1 = randperm(n_values);
%                     p2 = randperm(n_values);
%                     a = p1 - p2;
%                 end
%                 D(i1:i2) = all_ef_vect(c1,p1) - all_ef_vect(c2,p2); % differences
%                 i1 = i2 + 1;
%                 i2 = i2 + n_values;
%             end
%         end
%     end
% end




nbins = 100;
% histfit(D,nbins,'kernel')
% n = hist(D,nbins);
% maxn = max(n);
pd_D = fitdist(D,'kernel');
xpd = min(D):1/nbins:max(D);
ypd = pdf(pd_D,xpd);
csypd = cumsum(ypd);
csypd = csypd ./ max(csypd);
csypd = (csypd .* 2) - 1;

fprintf(' done\n')

%% Compute general preference map & prob map
Prefs = zeros(dims);
Probs = zeros(dims);
for i = 1:size(Prefs,1)
    for j = 1:size(Prefs,2)
        for k = 1:size(Prefs,3)
            if all_ef(i,j,k,1) ~= 0
                [a,ind] = sort(all_ef(i,j,k,:),4,'descend');
                a = squeeze(a);
                ind = squeeze(ind);
                Prefs(i,j,k) = ind(1);
                diff_SW = a(1) - a(2);
                [~,I] = min(abs(xpd - diff_SW));
                Probs(i,j,k) = csypd(I);
            end
        end
    end
end


