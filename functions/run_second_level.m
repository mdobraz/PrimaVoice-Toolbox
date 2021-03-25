ncons = numel(contrast_names);

if ~exist(paths.results_multi,'dir');mkdir(paths.results_multi);end % create folder if non-existant

%% Write log file
log_file_results = fullfile(paths.results_multi,['Results_second_level_analysis_' paths.analysis_name '.log']);
fid = fopen(log_file_results,'w');
fprintf(fid,'Subjects: ');
fprintf(fid,'%s ',paths.subjects{:});
fprintf(fid,'\n\n');


%% Combine runs fixed effects for each contrast %%%
MS_which_stats = '_t';
contrast = 1;
n_subj_to_combine = numel(paths.subjects);

p_val_peak = 0.05;
mags = {'t';'ef';'sd'};

for con = 1:ncons

    X = ones(n_subj_to_combine,1);

    infiles.t = cell(n_subj_to_combine,1);
    infiles.ef = cell(n_subj_to_combine,1);
    infiles.sd = cell(n_subj_to_combine,1);
    df_data = nan(1,n_subj_to_combine); % the row vector of degrees of freedom of the input files
    fwhm_data = nan(1,n_subj_to_combine);
    peak_thres = nan(1,n_subj_to_combine);
    clust_thres = nan(1,n_subj_to_combine);
    ext_thres = nan(1,n_subj_to_combine);
    tozip_in = zeros(n_subj_to_combine,numel(mags));
    for p = 1:n_subj_to_combine
        thres_file = fullfile(paths.results{p},['sub-' paths.subjects{p} '_res-' paths.results_name{p} '_' contrast_names{con} '_threshold.mat']);
        thres = load(thres_file);
        df_data(p) = thres.df.t;
        fwhm_data(p) = thres.fwhm_data;
        peak_thres(p) = thres.peak_threshold(thres.p_val_peak == p_val_peak);
        clust_thres(p) = thres.Cluster_threshold;
        ext_thres(p) = thres.extent_threshold_vol;
        for m = 1:numel(mags)
            infile = fullfile(paths.results{p},['In-' paths.template_name '_sub-' paths.subjects{p} '_res-' paths.results_name{p} '_' contrast_names{con} '_' mags{m} '.nii.gz']);
            [~,~,ext] = fileparts(infile);
            if strcmp(ext,'.gz')
                fprintf('Unzipping: %s\n',infile)
                system(['gunzip ' infile]);
                tozip_in(p,m) = 3;
            end
            infiles.(mags{m}){p} = infile(1:end-tozip_in(p,m));
        end
    end
    infiles_ef = char(infiles.ef);
    infiles_sd = char(infiles.sd);

    %% Brain volume for threshold
    [volume, number_voxels] = mask_vol(infiles.sd{1}); % "As input this function requires a mask_file which is the original Motion Corrected fMRI data"


    %% Fixed effects
    out_base = fullfile(paths.results_multi,['second_level_' paths.analysis_name '_fixed_' contrast_names{con}]);
    delete([out_base '*']);
    df = multistat(infiles_ef,infiles_sd,df_data,[],X,contrast,out_base,MS_which_stats,Inf); %tata!
    % Last argument of multistat, FWHM_VARATIO:
    %  - 0 will do no smoothing, and give a purely random effects analysis;
    %  - Inf will do complete smoothing to a global ratio of one, giving a purely fixed effects analysis.

    % find threshold
    [peak_threshold,extent_threshold_vol,Cluster_threshold] = stat_threshold(volume,number_voxels,mean(fwhm_data),df.t,p_val_peak);
    hd = fmris_read_image(infiles.sd{1});
    extent_threshold_vox = round(extent_threshold_vol / abs(prod(hd.vox)));
    save([out_base '_threshold.mat'],'peak_threshold','extent_threshold_vol','extent_threshold_vox','Cluster_threshold','df','p_val_peak')
     
    
    % Get Activations
    t_map_file = sprintf('%s_t.nii',out_base);
    [n_vox_mask,n_vox_maskP,meanT,meanTP,maxT] = get_activation(t_map_file,contrast_names{con},p_val_peak(1),p_val_peak,peak_threshold,extent_threshold_vox,Cluster_threshold,paths,0);
    
    fprintf(fid,'\n### AC activation for %s (fixed effects):\n',contrast_names{con});
    fprintf(fid,'Spatial extent (total) = %i\n',n_vox_mask);
    fprintf(fid,'Spatial extent (peaks) = %i\n',n_vox_maskP);
    fprintf(fid,'Mean t-value (total) = %g\n',meanT);
    fprintf(fid,'Mean t-value (peaks) = %g\n',meanTP);
    fprintf(fid,'Max t-value = %g\n',maxT);
     
    
    %% Random effects
    out_base = fullfile(paths.results_multi,['second_level_' paths.analysis_name '_random_' contrast_names{con}]);
    delete([out_base '*']);
    df = multistat(infiles_ef,infiles_sd,df_data,[],X,contrast,out_base,MS_which_stats,0); %tata!

    % save threshold
    save([out_base '_threshold.mat'],'peak_threshold','extent_threshold_vol','extent_threshold_vox','Cluster_threshold','df','p_val_peak')

    % find threshold
    % [peak_threshold,extent_threshold_vol,Cluster_threshold] = stat_threshold(volume,number_voxels,realign_flags.fwhm,df.t,p_val_peak);
    % extent_threshold_vox = round(extent_threshold_vol / abs(prod([scan_info(1,1).VoxelSizeX scan_info(1,1).VoxelSizeY scan_info(1,1).SliceThickness])));
    % save([out_base '_threshold.mat'],'peak_threshold','extent_threshold_vol','extent_threshold_vox','Cluster_threshold','df','p_val_peak')

    % Get Activations
    t_map_file = sprintf('%s_t.nii',out_base);
    [n_vox_mask,n_vox_maskP,meanT,meanTP,maxT] = get_activation(t_map_file,contrast_names{con},p_val_peak(1),p_val_peak,peak_threshold,extent_threshold_vox,Cluster_threshold,paths,0);
    
    fprintf(fid,'\n### AC activation for %s (random effects):\n',contrast_names{con});
    fprintf(fid,'Spatial extent (total) = %i\n',n_vox_mask);
    fprintf(fid,'Spatial extent (peaks) = %i\n',n_vox_maskP);
    fprintf(fid,'Mean t-value (total) = %g\n',meanT);
    fprintf(fid,'Mean t-value (peaks) = %g\n',meanTP);
    fprintf(fid,'Max t-value = %g\n',maxT);



    %% Coincidence map
    if ~stim_level
        out_base = fullfile(paths.results_multi,['second_level_' paths.analysis_name '_' contrast_names{con}]);
        delete([out_base '*']);
        coincidence_map(infiles.t,paths,peak_thres,ext_thres,clust_thres,out_base)

        out_base = fullfile(paths.results_multi,['second_level_' paths.analysis_name '_' contrast_names{con} '_inverse']);
        coincidence_map(infiles.t,paths,peak_thres,ext_thres,clust_thres,out_base,1)
    end

    %% Zip if needed
    for p = 1:n_subj_to_combine
        for m = 1:numel(mags)
            if tozip_in(p,m)
                fprintf('Zipping: %s\n',infiles.(mags{m}){p})
                system(['gzip ' infiles.(mags{m}){p}]);
            end
        end
    end
    
    

end

fprintf('\nZipping everything...')
system(sprintf('gzip %s/*.nii',paths.results_multi));
fprintf(' done.\n')


fprintf(fid,'\n### Thresholds:\n');
for i = 1:length(p_val_peak)
    fprintf(fid,'p = %g: Peak threshold = %g\n',p_val_peak(i),peak_threshold(i));
end
fprintf(fid,'Cluster threshold = %g\n',Cluster_threshold);
fprintf(fid,'Extent threshold = %g mm^3\n',extent_threshold_vol);
fprintf(fid,'Extent threshold = %g voxels (in the functional space)\n',extent_threshold_vox);


fclose(fid);
