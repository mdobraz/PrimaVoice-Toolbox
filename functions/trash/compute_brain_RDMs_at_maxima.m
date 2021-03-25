function compute_brain_RDMs_at_maxima(RDM_params,paths,ScanLog,coreg_params)
% Compute RDMs for each local maxima using voxels around the peak value.
%
% Parameters
% ----------
% mask_contrast: name of the contrast that will be used to mask the data
% peak_threshold: value at which this contrast t-map will thresholded
% mask_with_grey: if set to 1, non-grey matter voxels will be removed from analysis
% min_clust_size: only keep clusters from the contrast t-map that are larger than this value
% n_voxels: size of the max cluster. Examples:
            % 1: maximum only, 7: +faces, 19: +edges, 27: +corners
            % 33: +faces2, 57: +edges2, 81: +corners2, 125: 5-cube
% n_peaks: % nb of maxima to keep in each cluster.
           % if >= 1, same nb of max in each cluster
           % if < 1, the nb of peaks of each cluster will be:
           % floor((sizeR / n_voxels) * n_peaks), where sizeR is the 
           % size of the cluster being analyzed
% lastarg: either coreg_params or fixed_random (char, use fixed effects or random effects group analysis)


mask_contrast = RDM_params.contrast;
peak_threshold = RDM_params.peak_threshold;
mask_with_grey = RDM_params.mask_with_grey;
min_clust_size = RDM_params.min_clust_size;
n_voxels = RDM_params.n_voxels;
n_peaks = RDM_params.n_peaks;
if isfield(paths,'subjects')
    fixed_random = RDM_params.fixed_random;
    FRstr = [fixed_random '_'];
else
    FRstr = '';
end

if n_voxels > min_clust_size
    min_clust_size = n_voxels;
    warning('Needed neighbors is lower than minimal cluster extent. min_clust_size is now %d.\n', n_voxels);
end



if ischar(lastarg)
    fixed_random = lastarg;
elseif isstruct(lastarg)
    coreg_params = lastarg;
end
        

%% load t_file (or any file that will be used to find clusters and maxima)
if isfield(paths,'CV_results_sub_folder') && ~isempty(paths.CV_results_sub_folder)
    if isfield(paths,'subjects')
        path_t_maps = fullfile(paths.dataset,['second_level_analysis_' paths.CV_results_sub_folder]);
    else
        path_t_maps = fullfile(paths.results,paths.CV_results_sub_folder);
    end
    isCV = 1;
else
    isCV = 0;
    if isfield(paths,'subjects')
        path_t_maps = fullfile(paths.dataset,['second_level_analysis_' paths.analysis_name]);
    else
        path_t_maps = paths.results_multi;
    end
end


if isfield(paths,'subjects')
    if isCV
        t_file = fullfile(path_t_maps,['second_level_' paths.CV_results_sub_folder '_' fixed_random '_' mask_contrast '_t.nii.gz']);
    else
        t_file = fullfile(path_t_maps,['second_level_' paths.analysis_name '_' fixed_random '_' mask_contrast '_t.nii.gz']);
    end
else
    t_file = fullfile(path_t_maps,['sub-' paths.subject '_res-' paths.results_name '_' mask_contrast '_t.nii.gz']);
end
Pt = spm_vol(t_file);
Yt = spm_read_vols(Pt);

% if MION, inverse values
if strcmpi(paths.contrast_agent,'MION')
    Yt = -Yt;
end

%% create mask
mask = false(size(Yt));
mask(Yt >= peak_threshold) = true;


%% create mask
% mask = false(size(Yt));

% if ~extent_thres
%     % mask using the peak threshold
%     mask(Yt >= peak_threshold) = true;
% else
%     % mask using the cluster threshold
%     mask(Yt > Cluster_threshold) = true;
    
%     % remove clusters smaller than extent threshold in mask
%     [labels,num] = bwlabeln(mask,18); % 6-connected, 18-connected or 26-connected
%     for n = 1:num
%         cluster_size = sum(labels(:) == n);
%         if cluster_size < extent_threshold_vox
%             mask(labels == n) = false;
%         end
%     end
% end




%% mask with grey
if mask_with_grey
    % load tissue maps
    tissues = {'grey';'white';'csf'};
    for i = 1:numel(tissues)
        if isfield(paths,'subjects')
            Ptissues.(tissues{i}) = spm_vol(paths.template_tissue.(tissues{i})); % load from template space
        else
            Ptissues.(tissues{i}) = spm_vol(paths.func.tissues.(tissues{i})); % load from subject func space
        end
        Ytissues.(tissues{i}) = spm_read_vols(Ptissues.(tissues{i}));
    end

    % grey mask
    grey_mask = false(size(Yt));
    grey_mask(Ytissues.grey >= 0.1) = true; % add grey voxels
    grey_mask(Ytissues.white >= 0.1) = false; % remove white voxels
    grey_mask(Ytissues.csf >= 0.1) = false; % remove white voxels

    mask(~grey_mask) = false; % remove non-grey voxels
end


%% Segment
[tmp_labels,num] = bwlabeln(mask,18); % spm_max works with 26-connected regions
clust_sizes = nan(num,1);
for n = 1:num
    clust_sizes(n) = sum(tmp_labels(:) == n);
    % if cluster_size < min_clust_size
    %     mask(labels == n) = false;
    % end
end

[~,I] = sort(clust_sizes,'descend');

num = sum(clust_sizes >= min_clust_size); % remove clusters that have less than 3 voxels
clust_sizes = clust_sizes(I(1:num));
labels = zeros(size(tmp_labels));
for clust = 1:num
    labels(tmp_labels == I(clust)) = clust;
end
mask = logical(labels);


if ~num
    error('No Cluster found for these parameters.');
elseif num == 1
    fprintf('\n##### %i Cluster found in ''%s''\n',num,mask_contrast)
else
    fprintf('\n##### %i Clusters found in ''%s''\n',num,mask_contrast)
end


%% Get local maxima
% Identify where are values
indexes = find(mask == true);
[Lx, Ly, Lz] = ind2sub(size(mask),indexes);
vox_locs = [Lx'; Ly'; Lz'];

% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in voxels}
% A     - region number
[N, Z, M, A] = spm_max(Yt(indexes),vox_locs);

if ~length(A)
    error('No Maximum found for these parameters.');
elseif length(A) == 1
    fprintf('### %i Maximum in total\n',length(A))
else
    fprintf('### %i Maxima in total\n',length(A))
end

%% Create patch of n_voxels around each maximum
nmax = length(Z);
L = zeros(nmax,n_voxels);
for p = 1:nmax
    % Select all voxels included in the cluster of the peak
    sel_mask = zeros(size(Yt)); 
    sel_mask(labels == labels(M(1,p),M(2,p),M(3,p))) = 1;
    sel_indexes = find(sel_mask > 0);
    nvox = length(sel_indexes);

    % Convert voxel indexes into coordinates
    [vx_x, vx_y, vx_z] = ind2sub(size(Yt), sel_indexes);
    sel_locs = [vx_x, vx_y, vx_z];
    
    % Compute distance between each voxel in the peak
    sel_locs = sel_locs - repmat(M(:,p), 1, nvox)';
    sel_dists = zeros(nvox, 1);
    for v = 1:nvox
        sel_dists(v) = norm(sel_locs(v, :));
    end
    
    % Sort by distance and store only first voxels 
    [~,sort_idx] = sort(sel_dists);
    sorted_indexes = sel_indexes(sort_idx);
    L(p, :) = sorted_indexes(1:n_voxels);
end




%% Remove patches that have common voxels (keep the highest ones)
cond_triang = tril(~diag(1:size(L,1)));
pairs = nan(sum(cond_triang(:)),2);
[pairs(:,1),pairs(:,2)] = ind2sub(size(cond_triang),find(cond_triang));


torm = [];
for i = 1:size(pairs,1)
    p = pairs(i,1);
    q = pairs(i,2);
    if any(ismember(L(p,:),L(q,:)))
        id = [p q];
        [~,im] = min([Z(p) Z(q)]);
        % torm(i) = id(im);
        if ~ismember(id(im),torm)
            torm = [torm id(im)];
        end
    end
end
if ~isempty(torm)
    L(torm,:) = [];
    N(torm) = [];
    Z(torm) = [];
    A(torm) = [];
end

if ~length(A)
    error('No Maximum found for these parameters.');
elseif length(A) == 1
    fprintf('### %i Maximum after removing overlapping maxima\n',length(A))
else
    fprintf('### %i Maxima after removing overlapping maxima\n',length(A))
end




%% Filter to keep only n_peaks peaks by cluster
if n_peaks ~= 0 % if n_peaks == 0 take all the available maxima
    for a = 1:num
        if n_peaks < 1
            sizeR = unique(N(A==a));
            n_peaksR = floor((sizeR / n_voxels) * n_peaks);
            if ~n_peaksR; n_peaksR = 1; end
        else
            n_peaksR = n_peaks;
        end
        n_peaksR = min([n_peaksR;sum(A==a)]);
        if sum(A==a) > n_peaksR % remove non selected maxima
            p_idx = find(A==a);
            [~,I] = sort(Z(A==a),'descend');
            Z(p_idx(I(n_peaksR+1:end))) = [];
            A(p_idx(I(n_peaksR+1:end))) = [];
            N(p_idx(I(n_peaksR+1:end))) = [];
            L(p_idx(I(n_peaksR+1:end)),:) = [];
        end
    end
end

if ~length(A)
    error('No Maximum found for these parameters.');
elseif length(A) == 1
    fprintf('### %i Maximum after aiming for %g max in each cluster\n\n',length(A),n_peaks)
else
    fprintf('### %i Maxima after aiming for %g max in each cluster\n\n',length(A),n_peaks)
end

nclust = length(A);



%% sort L by Z
[~,I] = sort(Z,'descend');
L = L(I,:);

  
%% sorted labels
sorted_labels = zeros(size(Yt));
for clust = 1:nclust
    sorted_labels(L(clust,:)) = clust;
end


%% paths
thr_str = sprintf('%g',peak_threshold);
thr_str = strrep(thr_str,'.','p');
peaks_str = sprintf('%g',n_peaks);
peaks_str = strrep(peaks_str,'.','p');

if mask_with_grey
    grey_str = '_grey-masked';
else
    grey_str = '';
end
if isCV
    CV_str = ['_CV-' paths.CV_results_sub_folder];
else
    CV_str = '';
end

if ~isfield(paths,'L') || isempty(paths.L)
    L_str = sprintf('_L%i',numel(fieldnames(ScanLog.fMRISTAT)));
else
    L_str = sprintf('_L%i',paths.L);
end

RDM_name = sprintf('%sbrain_RDM_%s_thres-%s_clust-size-%i_n-vox-%i_n-peaks-%s%s%s%s',FRstr,mask_contrast,thr_str,min_clust_size,n_voxels,peaks_str,grey_str,CV_str,L_str);

RDMs_dir = fullfile(paths.results_multi,'brain_RDMs');
if ~exist(RDMs_dir,'dir');mkdir(RDMs_dir);end % create tmp folder if non-existant
RDM_dir = fullfile(RDMs_dir,RDM_name);
if ~exist(RDM_dir,'dir');mkdir(RDM_dir);end % create tmp folder if non-existant
delete(fullfile(RDM_dir,'*'));

if isfield(paths,'subjects')
    prefix = ['second_level_' paths.analysis_name '_' RDM_name];
else
    prefix = ['sub-' paths.subject '_res-' paths.results_name '_' RDM_name];
end


%% save labels volume
P = Pt;
label_names{1} = [prefix '_cluster-labels.nii'];
P.fname = fullfile(RDM_dir,label_names{1});
P.descrip = 'Labels';
P = rmfield(P,'pinfo');
spm_write_vol(P,labels);
system(['gzip ' P.fname]);

%% save sorted_labels volume
label_names{2} = [prefix '_maxima-labels.nii'];
P.fname = fullfile(RDM_dir,label_names{2});
spm_write_vol(P,sorted_labels);
system(['gzip ' P.fname]);



%% Register to the anat & template if subject analysis
if ~isfield(paths,'subjects') % subject analysis
    setenv('ANTSPATH',paths.ANTS_path);
    fprintf('Applying spatial transforms...\n\n')
    if strcmp(coreg_params.method,'bbr')
        t_average_to_anat_warp = '';
    elseif strcmp(coreg_params.method,'ants')
        t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
    else
        error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
    end

    for i = 1:2
        % to anat
        in_file = fullfile(RDM_dir,[label_names{i} '.gz']);
        out_file = fullfile(RDM_dir,['In-' paths.anat.name '_' label_names{i} '.gz']);
        system(sprintf('%s -i %s -r %s -o %s%s -t %s -n NearestNeighbor',...
            fullfile('$ANTSPATH','antsApplyTransforms'),...
            in_file,paths.anat.full,out_file,...
            t_average_to_anat_warp,paths.average_to_anat_xfm));
     
        % to template
        out_file = fullfile(RDM_dir,['In-' paths.template_name '_' label_names{i} '.gz']);
        system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n NearestNeighbor',...
            fullfile('$ANTSPATH','antsApplyTransforms'),...
            in_file,paths.anat.full,out_file,...
            paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
            t_average_to_anat_warp,paths.average_to_anat_xfm));
    end
end










%% cluster variable for cluster position
cluster = struct('maxT', cell(1,nclust),...
    'meanT', cell(1,nclust),...
    'max_coord', cell(1,nclust),...
    'Tvalues', cell(1,nclust),...
    'size', cell(1,nclust));


%% Get stim_names
fn_L = fieldnames(ScanLog.fMRISTAT);
if ~isfield(paths,'L') || isempty(paths.L)
    stim_names = ScanLog.fMRISTAT.(fn_L{end}).names(1:end-1); % do not take silence condition
else
    stim_names = ScanLog.fMRISTAT.(fn_L{paths.L}).names(1:end-1); % do not take silence condition
end
n_stims = length(stim_names);

% %% Create the pairs to compare
% cond_triang = tril(~diag(1:n_stims));
% pairs = nan(sum(cond_triang(:)),2);
% [pairs(:,1),pairs(:,2)] = ind2sub(size(cond_triang),find(cond_triang));



RDM_norm = zeros(n_stims,n_stims,nclust);

%% Load all stim maps
t_maps = cell(n_stims,1);

for i = 1:n_stims
    if isfield(paths,'subjects')
        t_file = fullfile(paths.results_multi,['second_level_' paths.analysis_name '_' fixed_random '_' stim_names{i} '_vs_silence_t.nii.gz']);
    else
        t_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' stim_names{i} '_vs_silence_t.nii.gz']);
    end
    P = spm_vol(t_file);
    t_maps{i} = spm_read_vols(P);
    if strcmpi(paths.contrast_agent,'MION')
        t_maps{i} = -t_maps{i};
    end
end


%% for each cluster
for clust = 1:nclust
    fprintf('\n### Cluster %i / %i\n',clust,nclust)
    fprintf('Computing bootstrapped CIs:\n')

    cluster(clust).maxT = max(Yt(sorted_labels == clust));
    mi = find(Yt == cluster(clust).maxT);
    [a,b,c] = ind2sub(size(Yt),mi);
    cluster(clust).max_coord = [a b c];
    cluster(clust).size = sum(sorted_labels(sorted_labels == clust));
    cluster(clust).meanT = mean(Yt(sorted_labels == clust));
    cluster(clust).Tvalues = Yt(sorted_labels == clust);

    cluster(clust).stim = struct('name', cell(1,n_stims),...
                                'Tvalues', cell(1,n_stims),...
                                'medianT', cell(1,n_stims),...
                                'median_ci1', cell(1,n_stims),...
                                'median_ci2', cell(1,n_stims),...
                                'median_L', cell(1,n_stims),...
                                'median_U', cell(1,n_stims),...
                                'meanT', cell(1,n_stims),...
                                'mean_ci1', cell(1,n_stims),...
                                'mean_ci2', cell(1,n_stims),...
                                'mean_L', cell(1,n_stims),...
                                'mean_U', cell(1,n_stims));
    thisActivityPattern = nan(n_stims,n_voxels);
    for i = 1:n_stims
        if ~mod(i,round(n_stims / 10))
            fprintf('%i%% ',round((i / n_stims) * 100))
        end
        cluster(clust).stim(i).name = stim_names{i};
        c = t_maps{i};
        c(sorted_labels ~= clust) = nan;
        x = c(~isnan(c(:)));
        cluster(clust).stim(i).Tvalues = x;
        [xpd,csypd] = boot_pd(10000,@median,x);
        [cluster(clust).stim(i).median_ci1,cluster(clust).stim(i).median_ci2,~,cluster(clust).stim(i).median_L,cluster(clust).stim(i).median_U] = CI_values(xpd,csypd,0.95);
        cluster(clust).stim(i).medianT = nanmedian(c(:));
        [xpd,csypd] = boot_pd(10000,@mean,x);
        [cluster(clust).stim(i).mean_ci1,cluster(clust).stim(i).mean_ci2,~,cluster(clust).stim(i).mean_L,cluster(clust).stim(i).mean_U] = CI_values(xpd,csypd,0.95);
        cluster(clust).stim(i).meanT = nanmean(c(:));
        thisActivityPattern(i,:) = x';
    end
    fprintf('\nConstructing RDM...')
    RDM(:,:,clust) = squareform(pdist(thisActivityPattern));%,userOptions.distance));
    fprintf('\n')
    % %% for each pair
    % for p = 1:size(pairs,1)
        
    %     if ~mod(p,round(size(pairs,1) / 10))
    %         fprintf('%i%% ',round((p / size(pairs,1)) * 100))
    %     end 
        
    %     c1 = t_maps{pairs(p,1)};
    %     c1(sorted_labels ~= clust) = nan;
    %     c1 = reshape(c1,numel(c1),1);
    %     c1(isnan(c1)) = [];
        
    %     c2 = t_maps{pairs(p,2)};
    %     c2(sorted_labels ~= clust) = nan;
    %     c2 = reshape(c2,numel(c2),1);
    %     c2(isnan(c2)) = [];
            
    %     RDM(pairs(p,1),pairs(p,2),clust) = norm(c1 - c2);
    % end
    % RDM(:,:,clust) = RDM(:,:,clust) + RDM(:,:,clust)';
    % fprintf('\n\n')
end





%% save RDM file
RDM_file = fullfile(RDM_dir,[prefix '.mat']);
save(RDM_file,'RDM','cluster','-v7');



%% Trash

% % Filter to keep only n_peaks peaks by clusters
% % For each cluster
% highest_Lpeaks = [];
% peaks_value = [];
% for a = 1:num
%     tmp_Z = Z(A==a);
%     tmp_M = M(:, A==a);
%     [sorted_Z, idx_Z] = sort(tmp_Z,'descend');
%     tmp_M = tmp_M(:, idx_Z);
%     ordered_peaks = sub2ind(size(Yt), tmp_M(1, :), tmp_M(2, :), tmp_M(3, :));
%     if n_peaks < 1
%         sizeR = unique(N(A==a));
%         n_peaksR = floor((sizeR / n_voxels) * n_peaks);
%         if ~n_peaksR; n_peaksR = 1; end
%     else
%         n_peaksR = n_peaks;
%     end
%     n_peaksR = min([n_peaksR;length(ordered_peaks)]);
    
%     highest_Lpeaks = [highest_Lpeaks, ordered_peaks(1:n_peaksR)];
%     peaks_value = [peaks_value, Yt(ordered_peaks(1:n_peaksR))];
% end
% [~,peaks_order] = sort(peaks_value,'ascend'); % lowest value is lowest maximum
% Lpeaks = highest_Lpeaks(peaks_order);
% nclust = length(Lpeaks);

% %% Search for voxels around each maximum
% L = zeros(nclust, n_voxels);
% for p = 1:nclust
%     % Select all voxels included in the cluster of the peak
%     sel_mask = zeros(size(Yt)); 
%     sel_mask(labels == labels(Lpeaks(p))) = 1;
%     sel_indexes = find(sel_mask > 0);
%     nvox = length(sel_indexes);

%     % Convert voxel indexes into coordinates
%     [vx_x, vx_y, vx_z] = ind2sub(size(Yt), sel_indexes);
%     sel_locs = [vx_x, vx_y, vx_z];
    
%     % Compute distance between each voxel in the peak
%     [xmax, ymax, zmax] = ind2sub(size(Yt), Lpeaks(p)); % location of peak p
%     sel_locs = sel_locs - repmat([xmax, ymax, zmax]', 1, nvox)';
%     sel_dists = zeros(nvox, 1);
%     for v = 1:nvox
%         sel_dists(v) = norm(sel_locs(v, :));
%     end
    
%     % Sort by distance and store only first voxels 
%     [sorted_dists,sort_idx] = sort(sel_dists);
%     sorted_indexes = sel_indexes(sort_idx);
%     L(p, :) = sorted_indexes(1:n_voxels);
% end







