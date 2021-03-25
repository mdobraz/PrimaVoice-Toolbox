if coreg_params.use_daily_anat
    %% find daily anat
    fprintf('\nFinding daily anatomic file...')
    load(paths.reference_scan_infos)
    daily_anat_file = spm_BIDS(BIDS,'data','sub',paths.subject,'type',paths.daily_anat_suffix,'ses',sprintf('%02.0f',Master.session),'modality','anat');
    if numel(daily_anat_file) > 1 % if several anat runs, average them
        daily_anat_runs = spm_BIDS(BIDS,'runs','sub',paths.subject,'type',paths.daily_anat_suffix,'ses',sprintf('%02.0f',Master.session),'modality','anat');
        if numel(daily_anat_runs) > 1 % if several anat runs, average them
            daily_anat_file = cell(numel(daily_anat_runs),1);
            for r = 1:numel(daily_anat_runs)
                a = spm_BIDS(BIDS,'data','sub',paths.subject,'type',paths.daily_anat_suffix,'ses',sprintf('%02.0f',Master.session),'modality','anat','run',daily_anat_runs{r});
                daily_anat_file{r} = a{1};
            end
            [da_path,da_name,ext] = fileparts(daily_anat_file{1});
            if strcmp(ext,'.gz'); [~,da_name,ext] = fileparts(da_name); end
            strtofind = ['_run-' daily_anat_runs{1}];
            idx = strfind(da_name,strtofind);
            da_name(idx:idx+length(strtofind)-1) = '';
            merge_out_file = fullfile(da_path,[da_name ext]);
            
            if exist(merge_out_file,'file') == 2; delete(merge_out_file);end
            merge_in_files = sprintf('%s ',daily_anat_file{:});
            system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,merge_out_file,merge_in_files)); % merge anat runs
            system(sprintf('gunzip %s',merge_out_file));
            
            mean_out_file = fullfile(da_path,[da_name 'Mean' ext]);
            if exist(mean_out_file,'file') == 2; delete(mean_out_file);end
            system(sprintf('%sfslmaths %s -Tmean %s',paths.FSL_prefix,merge_out_file,mean_out_file)); % average merged anat runs
            system(sprintf('gunzip %s',mean_out_file));
            
            daily_anat_file = mean_out_file;
        else
            daily_anat_file = daily_anat_file{1};
        end
    elseif numel(daily_anat_file) == 0
        error('No anatomic file of type ''%s'' found in session %i for subject ''%s''.\nCheck the filename of your anatomic scan and/or modify ''paths.daily_anat_suffix'' in your parameters file',paths.daily_anat_suffix,Master.session,paths.subject)
    else
        daily_anat_file = daily_anat_file{1};
    end
    
    %% Move to the resliced folder
    [da_path,da_name,ext] = fileparts(daily_anat_file);
    if strcmp(ext,'.gz')
        system(sprintf('gunzip %s',daily_anat_file));
        daily_anat_file = fullfile(da_path,da_name);
        [~,da_name,ext] = fileparts(daily_anat_file);
    end
    daily_anat_file_in_data = daily_anat_file;
    daily_anat_file = fullfile(paths.resliced,[da_name ext]);
    copyfile(daily_anat_file_in_data,daily_anat_file);
    fprintf(' done.\n')
    
    %% preprocess daily anat (crop, N4, sanlm, BET)
    a = dir(fullfile(paths.resliced,[da_name '*_BET.nii'])); % check if a previous preprocessing already exists & ask whether the user wants to perform a new one or keep the old one
    perform_preproc = 1;
    if ~isempty(a)
        answer = questdlg(sprintf('A preprocessed daily anatomic file already exist: %s',a(1).name),'Replace old preprocessing?','Redo preprocessing','Keep old preprocessing','Redo preprocessing');
        if strcmp(answer,'Keep old preprocessing')
            perform_preproc = 0;
        end
    end
    
    if perform_preproc
        in_file = daily_anat_file; % in_file for strong_register_preprocessing script
        method = 'anat_manualBET'; % method for strong_register_preprocessing script
        strong_register_preprocessing
        unbet_da_file = in_file; % in_file preprocessed but unbetted
        BET_da_file = BET_out_file; % same but betted
    else % get filenames
        BET_da_file = fullfile(paths.resliced,a(1).name);
        unbet_da_file = BET_da_file;
        idx = strfind(unbet_da_file,'_BET');
        unbet_da_file(idx:idx+length('_BET')-1) = '';
    end
    
    %% FLIRT to main anat (using strong_register_FLIRT)
    flirt_in_file = BET_da_file; % the preprocessed BET daily anat file
    ref_file = paths.anat_file_brain;
    ref_brain_mask = paths.anat_file_brain_mask;
    [xfm_da2ma_file,xfm_ma2da_file,flirt_out_file,da_brain_file,da_brain_mask] = strong_register_FLIRT(unbet_da_file,flirt_in_file,ref_file,ref_brain_mask,coreg_params.daily_anat_flirt_cost,paths);
    
    % white transform preparation
    [da_brain_path,da_brain_name] = fileparts(da_brain_file);
    da_white_file = fullfile(da_brain_path,sprintf('%s_white.nii',da_brain_name));
    white_in_file = paths.anat_file_brain_segmentation.(paths.which_seg).tissues.white;
    if exist(da_white_file,'file') == 2; delete(da_white_file);end
    
    %% if FNIRT is wanted
    switch coreg_params.daily_anat_method
        case {'fnirt','fnirt_brains'}
            fprintf('Non-linear registration... (this may take several minutes...)')
            if strcmp(coreg_params.daily_anat_method,'fnirt')
                fnirt_in_file = unbet_da_file;
                fnirt_ref_file = paths.anat_file;
            elseif strcmp(coreg_params.daily_anat_method,'fnirt_brains')
                fnirt_in_file = da_brain_file;
                fnirt_ref_file = paths.anat_file_brain;
            end
            [fnirt_in_path,fnirt_in_name] = fileparts(fnirt_in_file);
            [~,fnirt_ref_name] = fileparts(fnirt_ref_file);
            fnirt_in_aff = xfm_da2ma_file; % the affine transform previously obtained with flirt
            fnirt_out_file = fullfile(fnirt_in_path,sprintf('%s_FNIRT_to_%s.nii',fnirt_in_name,fnirt_ref_name));
            fnirt_cout_file = fullfile(fnirt_in_path,sprintf('%s_FNIRT_to_%s_COEFFS.nii',fnirt_in_name,fnirt_ref_name)); % field coefficients
            system(sprintf('%sfnirt --ref=%s --in=%s --iout=%s --aff=%s --cout=%s',paths.FSL_prefix,fnirt_ref_file,fnirt_in_file,fnirt_out_file,fnirt_in_aff,fnirt_cout_file));
            fnirt_cout_inv_file = fullfile(fnirt_in_path,sprintf('%s_FNIRT_to_%s_COEFFS_inv.nii',fnirt_in_name,fnirt_ref_name)); % field coefficients
            system(sprintf('%sinvwarp -r %s -w %s -o %s',paths.FSL_prefix,fnirt_in_file,fnirt_cout_file,fnirt_cout_inv_file)); % inverse warp
            
            % apply warp to white matter
            system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,fnirt_in_file,white_in_file,da_white_file,fnirt_cout_inv_file));
        case 'flirt'
            % apply xfm to white matter
            system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,white_in_file,da_brain_file,da_white_file,xfm_ma2da_file));
        otherwise
            error('''%s'' method for daily anatomic registration does not exist. Available methods are:\n''flirt'', ''fnirt'', ''fnirt_brains''',coreg_params.daily_anat_method)
    end
    
    system(sprintf('gunzip %s',da_white_file));
end

%% Preprocess functional average
in_file = paths.average;
switch coreg_params.method
    case {'bbr_manual','flirt'}
        [~,av_name] = fileparts(paths.average);
        a = dir(fullfile(paths.resliced,[av_name '*_BET.nii'])); % check if a previous preprocessing already exists & ask whether the user wants to perform a new one or keep the old one
        perform_preproc = 1;
        if ~isempty(a)
            answer = questdlg(sprintf('A preprocessed Average file already exist: %s',a(1).name),'Replace old preprocessing?','Redo preprocessing','Keep old preprocessing','Redo preprocessing');
            if strcmp(answer,'Keep old preprocessing')
                perform_preproc = 0;
            end
        end
        
        if perform_preproc
            method = 'manualBET';
            strong_register_preprocessing
            unbet_av_file = in_file;
            BET_av_file = BET_out_file;
        else
            BET_av_file = fullfile(paths.resliced,a(1).name);
            unbet_av_file = BET_av_file;
            idx = strfind(unbet_av_file,'_BET');
            unbet_av_file(idx:idx+length('_BET')-1) = '';
        end
    case 'bbr'
        method = 'autoBET';
        strong_register_preprocessing
        BET_av_file = BET_out_file;
        unbet_av_file = paths.average;
    otherwise
        error('''%s'' method for coregistration does not exist. Available methods are:\n''bbr'', ''bbr_manual'', ''flirt''',coreg_params.method)
end

%% Coregister func to anat (either to the main anat or to the daily anat)
if coreg_params.use_daily_anat
    ref_file = unbet_da_file;
    ref_file_brain = da_brain_file;
    ref_brain_mask = da_brain_mask;
    wmseg_file = da_white_file;
else
    ref_file = paths.anat_file;
    ref_file_brain = paths.anat_file_brain;
    ref_brain_mask = paths.anat_file_brain_mask;
    wmseg_file = paths.anat_file_brain_segmentation.(paths.which_seg).tissues.white;
end

% [xfm_i2r_file,~,av_to_anat] = strong_register_FLIRT(unbet_av_file,BET_av_file,ref_file_brain,ref_brain_mask,coreg_params.flirt_cost,paths);

switch coreg_params.method
    case {'bbr','bbr_manual'} % using epi_reg (BBR method)
        [~,in_name,ext] = fileparts(unbet_av_file);
        if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end
        [~,ref_name,ext] = fileparts(ref_file);
        if strcmp(ext,'.gz'); [~,ref_name] = fileparts(ref_name); end
        
        av_to_anat = fullfile(paths.resliced,sprintf('%s_BBR_to_%s.nii',in_name,ref_name));
        if exist(av_to_anat,'file') == 2; delete(av_to_anat);end
        system(sprintf('%sepi_reg --wmseg=%s --epi=%s --t1=%s --t1brain=%s --out=%s',paths.FSL_prefix,wmseg_file,BET_av_file,ref_file,ref_file_brain,av_to_anat));
        %         xfm_i2r_bbr_file = fullfile(paths.resliced,sprintf('%s_BBR_to_%s.xfm',in_name,ref_name));
        %         system(sprintf('%sflirt -in %s -ref %s -omat %s -out %s -dof 6 -cost bbr -searchcost bbr -init %s -wmseg %s',paths.FSL_prefix,BET_av_file,ref_file_brain,xfm_i2r_bbr_file,av_to_anat,xfm_i2r_file,wmseg_file));
        %         system(sprintf('gunzip %s',av_to_anat));
        %         xfm_i2r_file = xfm_i2r_bbr_file;
        xfm_i2r_file = fullfile(paths.resliced,sprintf('%s_BBR_to_%s.xfm',in_name,ref_name));
        movefile(fullfile(paths.resliced,sprintf('%s_BBR_to_%s.mat',in_name,ref_name)),xfm_i2r_file)
        system(sprintf('gunzip %s',av_to_anat));
    case 'flirt' % using strong_register_FLIRT
        [xfm_i2r_file,~,av_to_anat] = strong_register_FLIRT(unbet_av_file,BET_av_file,ref_file_brain,ref_brain_mask,coreg_params.flirt_cost,paths);
end



if exist(paths.average_to_anat,'file') == 2; delete(paths.average_to_anat);end

%% Concatenate transforms / warps if daily anat was used
if coreg_params.use_daily_anat
    switch coreg_params.daily_anat_method
        case {'fnirt','fnirt_brains'}
            system(sprintf('%sconvertwarp --ref=%s --warp1=%s --premat=%s --out=%s --relout',paths.FSL_prefix,paths.anat_file,fnirt_cout_file,xfm_i2r_file,paths.average_to_anat_warp));
            system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s',paths.FSL_prefix,paths.anat_file,BET_av_file,paths.average_to_anat,paths.average_to_anat_warp)); % apply to average
            system(sprintf('%sinvwarp -r %s -w %s -o %s',paths.FSL_prefix,BET_av_file,paths.average_to_anat_warp,paths.anat_to_average_warp)); % inverse warp
        case 'flirt'
            system(sprintf('%sconvert_xfm -omat %s -concat %s %s',paths.FSL_prefix,paths.average_to_anat_xfm,xfm_da2ma_file,xfm_i2r_file));
            system(sprintf('%sflirt -in %s -ref %s -out %s -init %s -applyxfm',paths.FSL_prefix,BET_av_file,paths.anat_file,paths.average_to_anat,paths.average_to_anat_xfm)); % apply to average
    end
    system(sprintf('gunzip %s',paths.average_to_anat));
else
    movefile(xfm_i2r_file,paths.average_to_anat_xfm); % final transform is xfm_i2r_file
    movefile(av_to_anat,paths.average_to_anat);
end

system(sprintf('%sconvert_xfm -omat %s -inverse %s',paths.FSL_prefix,paths.anat_to_average_xfm,paths.average_to_anat_xfm)); % inverse xfm


%% Save figure of anat overlaid with template
view_slice_overlay(paths.anat_file,paths.average_to_anat,0)
figurewrite(fullfile(paths.resliced,'Average_to_anat')) % Using GLMdenoise function
view_slice_overlay(paths.average_to_anat,paths.anat_file_brain_segmentation.(paths.which_seg).tissues.white,0,[],0.3,[],lines(2))
figurewrite(fullfile(paths.resliced,'White_on_average_to_anat')) % Using GLMdenoise function




%% Apply xfm or warp on each tissue in order to put them in the average functional space
tissues = fieldnames(paths.anat_file_brain_segmentation.(paths.which_seg).tissues);
for i = 1:numel(tissues)
    in_file = paths.anat_file_brain_segmentation.(paths.which_seg).tissues.(tissues{i});
    out_file = paths.func_tissues.(tissues{i});
    if exist(out_file,'file') == 2; delete(out_file);end
    if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
        system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,paths.average,in_file,out_file,paths.anat_to_average_warp));
    else
        system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,in_file,paths.average,out_file,paths.anat_to_average_xfm));
    end
    system(sprintf('gunzip %s',out_file));
    
    probs = paths.tissues_probs;
    for p = 1:numel(probs)
        prob = probs{p};
        in_file = paths.(sprintf('anat_file_tissues_%sp',probs{p})).(tissues{i});
        out_file = paths.(sprintf('func_tissues_%sp',probs{p})).(tissues{i});
        if exist(out_file,'file') == 2; delete(out_file);end
        if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
            system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,paths.average,in_file,out_file,paths.anat_to_average_warp));
        else
            system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,in_file,paths.average,out_file,paths.anat_to_average_xfm));
        end
        system(sprintf('gunzip %s',out_file));
    end
end

%% Apply xfm on white + csf mask
in_file = paths.anat_file_tissues_white_csf_mask;
out_file = paths.func_tissues_white_csf_mask;
if exist(out_file,'file') == 2; delete(out_file);end
if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
    system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,paths.average,in_file,out_file,paths.anat_to_average_warp));
else
    system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,in_file,paths.average,out_file,paths.anat_to_average_xfm));
end
system(sprintf('gunzip %s',out_file));


view_slice_overlay(paths.average,out_file,0,[],0.3,[],lines(2))
figurewrite(fullfile(paths.resliced,'White_CSF_on_average')) % Using GLMdenoise function


%% Apply xfm on brain mask & STG-STS mask
in_file = paths.anat_file_brain_mask;
out_file = paths.func_brain_mask;
if exist(out_file,'file') == 2; delete(out_file);end
if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
    system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,paths.average,in_file,out_file,paths.anat_to_average_warp));
else
    system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,in_file,paths.average,out_file,paths.anat_to_average_xfm));
end
system(sprintf('gunzip %s',out_file));

view_slice_overlay(paths.average,out_file,0,[],0.3,[],lines(2))
figurewrite(fullfile(paths.resliced,'Average_brain')) % Using GLMdenoise function

if isfield(paths,'STS_STG_mask')
    in_file = paths.STS_STG_reg_mask;
    out_file = paths.func_STS_STG_reg_mask;
    if exist(out_file,'file') == 2; delete(out_file);end
    if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
        system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s --interp=nn',paths.FSL_prefix,paths.average,in_file,out_file,paths.anat_to_average_warp));
    else
        system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,in_file,paths.average,out_file,paths.anat_to_average_xfm));
    end
    system(sprintf('gunzip %s',out_file));
end

%% Create full filled volume (mask used when no mask is specified)
if exist(paths.no_mask,'file') == 2; delete(paths.no_mask);end
system(sprintf('%sfslmaths %s -thr 0 %s -odt input',paths.FSL_prefix,paths.average,paths.no_mask));
system(sprintf('gunzip %s',paths.no_mask));





