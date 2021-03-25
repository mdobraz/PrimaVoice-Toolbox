function compute_clustered_betas(p_val,conds,pref_cons,contrasts,dims,paths,coreg_params,extent_thres,name_suffix)

p_val_str = sprintf('%g',p_val);
p_val_str = strrep(p_val_str,'.','p');

if extent_thres
    thres_suffix = [p_val_str '_clust-thres'];
else
    thres_suffix = [p_val_str '_peak-thres'];
end

if ~exist('name_suffix','var')
    name_suffix = thres_suffix;
else
    name_suffix = [thres_suffix name_suffix];
end


%% Get betas from vs_silence contrasts
nconds = length(conds);
all_ef = nan([dims nconds]);
% all_ef_vect_cell = cell(nconds,1);
for i = 1:nconds
    % load ef_file
    ef_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{conds(i)} '_ef.nii.gz']);
    Pef = spm_vol(ef_file);
    Yef = spm_read_vols(Pef); % put all contrasts in one variable
    % if MION, inverse values
    if strcmpi(paths.contrast_agent,'MION')
        Yef = -Yef;
    end
    all_ef(:,:,:,i) = Yef;
%     all_ef_vect_cell{i} = Yef(logical(Yef));
end


%% Create files for clustered beta-values plots (display with plot_clustered_betas)
fprintf('Computing clustered beta-values for %s:\n',name_suffix)
iterations = 10000;

for icon = 1:length(pref_cons)
    con = pref_cons(icon);
    fprintf('\t%s...\n',contrasts.names{con})
    
    % load thres_file
    thres_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_threshold.mat']);
    load(thres_file)
    
    % load t-file
    t_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_t.nii.gz']);
    Pt = spm_vol(t_file);
    Yt = spm_read_vols(Pt);
    
    % load ef_file
    ef_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_ef.nii.gz']); 
    Pef = spm_vol(ef_file);
    Yef = spm_read_vols(Pef);
    
    % if MION, inverse values
    if strcmpi(paths.contrast_agent,'MION')
        Yef = -Yef;
        Yt = -Yt;
    end
    
    % create mask
    mask = false(size(Yt));

    if ~extent_thres
        % mask using the peak threshold
        mask(Yt >= peak_threshold(p_val_peak == p_val)) = true;
    else
        % mask using the cluster threshold
        mask(Yt > Cluster_threshold) = true;

        % remove clusters smaller than extent threshold in mask
        [labels,num] = bwlabeln(mask,18); % 6-connected, 18-connected or 26-connected
        for n = 1:num
            cluster_size = sum(labels(:) == n);
            if cluster_size < extent_threshold_vox
                mask(labels == n) = false;
            end
        end
    end
    
    % label mask
    [labels,nclust] = bwlabeln(mask,18); % segmentation
 
    if nclust > 0
        % sort nums by size
        clust_sizes = nan(nclust,1);
        for clust = 1:nclust
            clust_sizes(clust) = sum(labels(:) == clust);
        end
        [~,I] = sort(clust_sizes,'descend');
        
        nclust = sum(clust_sizes > 2); % remove clusters that have less than 3 voxels
        
        sorted_labels = zeros(size(labels));
        for clust = 1:nclust
            sorted_labels(labels == I(clust)) = clust;
        end

        % get clusters betas
        Y_above_thresh = Yef;
        Y_above_thresh(~mask) = 0;
  
        cluster.maxT = nan(nclust,1);
        cluster.max_coord = cell(nclust,1);
        cluster.xpd = cell(nclust,length(conds));
        cluster.csypd = cell(nclust,length(conds));
        cluster.ypd = cell(nclust,length(conds));
        cluster.values = cell(nclust,length(conds));
        for clust = 1:nclust
            cluster.maxT(clust) = max(Yt(sorted_labels == clust));
            mi = find(Y_above_thresh == cluster.maxT(clust));
            [a,b,c] = ind2sub(size(Y_above_thresh),mi);
            cluster.max_coord{clust} = [a b c];
            
            clust_std_ef = all_ef;
            for i = 1:length(conds)
                tmp_ef = squeeze(all_ef(:,:,:,i));
                tmp_ef(sorted_labels ~= clust) = nan; % remove voxels outside of cluster
                clust_std_ef(:,:,:,i) = tmp_ef;
                cluster.values{clust,i} = tmp_ef(~isnan(tmp_ef(:)));
            end
            medianC = median(clust_std_ef(~isnan(clust_std_ef(:))));
            stdC = std(clust_std_ef(~isnan(clust_std_ef(:))));
            
            for i = 1:length(conds) % compute pd
                tmpC = cluster.values{clust,i};
%                 tmpC = (tmpC - medianC) ./ stdC;
                [cluster.xpd{clust,i},cluster.csypd{clust,i},cluster.ypd{clust,i}] = boot_pd(iterations,@median,tmpC);
            end
        end
        % save cluster file
        clust_betas_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_clusters_betas_' name_suffix '.mat']);
        save(clust_betas_file,'cluster');
        
        % save sorted labels volume
        P = Pef;
        P.fname = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_clusters_betas_' name_suffix '_labels.nii']);
        P.descrip = 'Labels';
        P = rmfield(P,'pinfo');
        spm_write_vol(P,sorted_labels);
        
        % Register it to the anat
        [~,out_name] = fileparts(P.fname);
        out_file = fullfile(paths.results_multi,['Registered_' out_name '.nii']);
        if exist(out_file,'file') == 2; delete(out_file);end
        if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
            system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp nn',paths.FSL_prefix,paths.anat_file,P.fname,out_file,paths.average_to_anat_warp));
        else
            system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour',paths.FSL_prefix,P.fname,paths.anat_file,out_file,paths.average_to_anat_xfm));
        end
    end
end
