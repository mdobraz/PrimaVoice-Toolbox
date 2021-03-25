%% SCRIPT
% to be run after parameters_*

%% Run spm_BIDS
fprintf('Retrieving BIDS informations...')
BIDS = spm_BIDS(paths.dataset);
fprintf(' done.\n\n')

%% Get species
known_species = {'human';'macaque';'marmoset';'baboon'};

if isempty(BIDS.participants)
    error('There must be a ''participants.tsv'' file in the BIDS folder!!')
elseif ~isfield(BIDS.participants,'species')
    error('There must be a ''species'' column in the ''participants.tsv'' file!!')
else
    pn = find(strcmp(BIDS.participants.participant_id,['sub-' paths.subject]));
    if isempty(pn)
        error('Could not find participant in the ''participants.tsv'' file!!')
    elseif length(pn) > 1
        error('There are more than one participant with the ID ''%s'' in the ''participants.tsv'' file!!',['sub-' paths.subject])
    else
        paths.species = BIDS.participants.species{pn};
        if ~strcmp(paths.species,known_species)
            error('Species in the ''participants.tsv'' file must be exactly one of those:\n%s',sprintf('%s\n',known_species{:}))
        end
        paths.contrast_agent = BIDS.participants.contrast_agent{pn};
        
        fprintf('Subject: %s\n',BIDS.participants.name{pn})
        if ~strcmpi(paths.contrast_agent,'none')
            fprintf('Contrast agent: %s\n',BIDS.participants.contrast_agent{pn})
        end
    end
end

%% Unzip files if necessary
BIDS = gunzip_all_data(BIDS,paths,1);


%% Paths

if ~isfield(paths,'analysis_folder')
    paths.analysis_folder = paths.dataset;
    fprintf('\nAnalysis will be done in the dataset folder\nAdd a''paths.analysis_folder'' to specify another analysis folder\n')
end

paths.subject_analysis = fullfile(paths.analysis_folder,['analysis_sub-' paths.subject]); % All analyses for this subject

paths.anat_dir = fullfile(paths.subject_analysis,'anat'); % folder where the segmented brain files will be

paths.segmentation_root = fullfile(paths.anat_dir,'segmentation'); % folder where the segmented brain files will be

paths.analysis = fullfile(paths.subject_analysis,paths.realign_name); % folder where one analysis corresponding to a certain realign name will be

paths.ME = sprintf('%s/T2star_maps',paths.subject_analysis); % Multi echo results

paths.realign = sprintf('%s/realignment',paths.analysis); % realignment matrices

paths.tmp_nii = sprintf('%s/tmp_nii',paths.analysis); % tmp folder to store original bols with updated headers

paths.rp_files = sprintf('%s/rp_files',paths.realign); % folder for the rp_files

paths.resliced = sprintf('%s/resliced',paths.analysis); % folder where the resliced motion corrected images will be

paths.resliced_fig = sprintf('%s/figures',paths.resliced); % folder for the coreg figures

paths.fieldmaps = sprintf('%s/fieldmaps',paths.analysis); % folder where the all the fieldmaps associated images will be

paths.PCA = fullfile(paths.analysis,'PCA'); % folder where the PCA results will be

paths.T1_normalized = fullfile(paths.analysis,'T1_normalized'); % folder where the PCA results will be

paths.results = fullfile(paths.analysis,['results_' paths.results_name]); % folder where the results will be

paths.GLMdenoise = fullfile(paths.results,sprintf('GLMdenoise_%i-white-PCs_%i-csf-PCs',GLM_params.GLMdenoise_nPCs.white,GLM_params.GLMdenoise_nPCs.csf)); % folder where the PCA results will be

paths.AC_folder = fullfile(paths.results,'Activations');

paths.costmatrix = fullfile(paths.resliced,'cost_matrix.mat');

paths.bash = fullfile(fileparts(which('RunScript.m')),'functions','bash');


%% create sessions / runs structure
SR = struct_sess_run(BIDS,sessions,runs,paths.subject,paths.task);

paths.in_jack = 0;
ST = dbstack;
for i = 1:numel(ST)
    if strcmp(ST(i).name,'launch_jackknife')
        paths.in_jack = 1;
        [sessions,runs,rem_session,rem_run] = get_sessions_runs_jackknife(SR,run_index);
        SR = struct_sess_run(BIDS,sessions,runs,paths.subject,paths.task);
        paths.results_jack_name = sprintf('%s_Jackknife_ses-%02.0f-run-%02.0f',paths.results_name,rem_session,rem_run);
        paths.results_multi = fullfile(paths.results,sprintf('Jackknife_%02.0f-%02.0f',rem_session,rem_run));
        break
    end
    if strcmp(ST(i).name,'select_jackrank')
        paths.in_jack = 1;
        jackknife_results
        [sessions,runs] = select_from_jack_ranks(sessions,runs,jack_ranks,max_rank);
        SR = struct_sess_run(BIDS,sessions,runs,paths.subject,paths.task);
        paths.results_jack_name = sprintf('%s_Jackknife_rank-%03.0f',paths.results_name,max_rank);
        paths.results_multi = fullfile(paths.results,sprintf('Jackknife_rank-%03.0f',max_rank));
        break
    end
end

if isempty(SR)
    fprintf('\n\n')
    warning('No functional data was found!!')
else
    n_total_runs = 0;
    for s = 1:numel(SR)
        n_total_runs = n_total_runs + length(SR(s).runs);
    end
end

%% End of Paths

if ~paths.in_jack
    if isfield(paths,'results_sub_folder') && ~isempty(paths.results_sub_folder)
        paths.results_multi = fullfile(paths.results,paths.results_sub_folder);
    else
        paths.results_multi = paths.results;
    end
end


%% Templates
templates % run templates.m


%% Main anat file - T1w
T1w_files_json = get_anat(BIDS,paths,paths.main_anat_suffix,paths.main_anat_session);
T1w_files = spm_jsonread(T1w_files_json);
anat_name = T1w_files.name;

%% T2 file
if isfield(paths,'T2_anat_suffix') && ~isempty(paths.T2_anat_suffix)
    T2w_files_json = get_anat(BIDS,paths,paths.T2_anat_suffix,paths.T2_anat_session);
    T2w_files = spm_jsonread(T2w_files_json);
else
    T2w_files = T1w_files;
end



%% Segmentation paths
paths.segmentation = fullfile(paths.segmentation_root,[paths.template_name '_' anat_name]); % folder where the segmented brain files will be
paths.preprocessing = fullfile(paths.segmentation,'preprocessing'); % folder where the preprocessing of the anat will happen
paths.seg_fig = fullfile(paths.segmentation,'figures'); % folder where the control figures will be stored
paths.surfaces = fullfile(paths.segmentation,'surfaces'); % folder where the surfaces will be stored



%% Segmentation files
tissues = {'grey';'white';'csf'};
paths.anat.tissue_probs = {'01';'50';'99'};
if isfield(T1w_files,paths.template_name) && isfield(T1w_files.(paths.template_name),'brain_segmented')
    paths.anat.full = T1w_files.(paths.template_name).fullhead;
    paths.anat.brain = T1w_files.(paths.template_name).brain;
    paths.anat.brain_mask = T1w_files.(paths.template_name).brain_mask;
    paths.anat.segmentation = T1w_files.(paths.template_name).brain_segmented;
    paths.anat.name = T1w_files.(paths.template_name).name;

    % tissues
    probs = paths.anat.tissue_probs;
    for i = 1:numel(tissues)
        paths.anat.tissues.(tissues{i}) = T1w_files.(paths.template_name).tissue_probs.(tissues{i});
        [~,anat_name] = fileparts(fileparts2(paths.anat.tissues.(tissues{i})));
        paths.func.tissues.(tissues{i}) = fullfile(paths.resliced,[anat_name '_in-func.nii']);
        for p = 1:numel(probs)
            paths.anat.tissues.(sprintf('p%s',probs{p})).(tissues{i}) = T1w_files.(paths.template_name).tissue_bins.(tissues{i}).(sprintf('p%s',probs{p}));
            [~,anat_name] = fileparts(fileparts2(paths.anat.tissues.(sprintf('p%s',probs{p})).(tissues{i})));
            paths.func.tissues.(sprintf('p%s',probs{p})).(tissues{i}) = fullfile(paths.resliced,[anat_name '_in-func.nii']);
        end
    end

    % Masks
    pca_tissues = {'white';'csf'};
    for i = 1:numel(pca_tissues)
        paths.anat.pca_mask.(pca_tissues{i}) = T1w_files.(paths.template_name).pca_masks.(pca_tissues{i});
        [~,anat_name] = fileparts(fileparts2(paths.anat.pca_mask.(pca_tissues{i})));
        paths.func.pca_mask.(pca_tissues{i}) = fullfile(paths.resliced,[anat_name '_in-func.nii']);
    end
    [~,anat_name] = fileparts(fileparts2(paths.anat.brain_mask));
    paths.func.brain_mask = fullfile(paths.resliced,[anat_name '_in-func.nii']);

    paths.no_mask = fullfile(paths.resliced,['sub-' paths.subject '_func_no_mask.nii']);
    mask_types = {'';'brain';'grey'};
    if ~any(strcmp(GLM_params.mask_type,mask_types))
        error('Mask type ''%s'' unknown',GLM_params.mask_type)
    end

    switch GLM_params.mask_type
        case ''
            paths.func_mask = paths.no_mask;
        case 'brain'
            paths.func_mask = paths.func.brain_mask;
            paths.anat_mask = paths.anat.brain_mask;
        case 'grey'
            paths.func_mask = paths.func.tissues.(sprintf('p%s',probs{1})).grey;
            paths.anat_mask = paths.anat.tissues.(sprintf('p%s',probs{1})).grey;
    end
else
    warning('No segmentation of the anatomic scan yet. Launch ''brain_segmentation'' and re-run parameters')
end
% spm_jsonread(get_corresp_json(anat_files.fullhead))


%% Template transforms
paths.ANTs_out_base = fullfile(paths.segmentation,'SyN_template_to_anat');
paths.temp_to_anat_xfm = [paths.ANTs_out_base '0GenericAffine.mat'];
paths.anat_to_temp_xfm = [paths.ANTs_out_base '0InverseGenericAffine.mat'];
paths.temp_to_anat_warp = [paths.ANTs_out_base '1Warp.nii.gz'];
paths.anat_to_temp_warp = [paths.ANTs_out_base '1InverseWarp.nii.gz'];



%% Func to anat registration files
reg_suffix = [paths.main_anat_suffix '_in-' paths.template_name];
paths.average_to_anat = fullfile(paths.resliced,sprintf('Average_%s_to_%s_using_%s.nii.gz',paths.realign_name,reg_suffix,coreg_params.method));
paths.average_to_anat_xfm = fullfile(paths.resliced,sprintf('Average_%s_to_%s_using_%s_xfm.mat',paths.realign_name,reg_suffix,coreg_params.method));
paths.anat_to_average_xfm = fullfile(paths.resliced,sprintf('%s_to_average_%s_using_%s_xfm.mat',reg_suffix,paths.realign_name,coreg_params.method));
paths.average_to_anat_warp = fullfile(paths.resliced,sprintf('Average_%s_to_%s_using_%s_warp.nii.gz',paths.realign_name,reg_suffix,coreg_params.method));
paths.anat_to_average_warp = fullfile(paths.resliced,sprintf('%s_to_average_%s_using_%s_warp.nii.gz',reg_suffix,paths.realign_name,coreg_params.method));

paths.average_to_temp = fullfile(paths.resliced,sprintf('Average_%s_to_%s_using_%s.nii.gz',paths.realign_name,paths.template_name,coreg_params.method));

%% Other files
paths.reference_scan = fullfile(paths.resliced,'Reference_scan.nii');
paths.reference_scan_infos = fullfile(paths.realign,'Reference_scan_infos.mat');
if paths.in_jack
    paths.log_file_params = sprintf('%s/Params_%s.log',paths.analysis,paths.results_jack_name);
else
    paths.log_file_params = sprintf('%s/Params_%s.log',paths.analysis,paths.results_name);
end
paths.log_file_realign = sprintf('%s/Realignment.log',paths.analysis);

paths.average = fullfile(paths.resliced,sprintf('Average_%s.nii',paths.realign_name));





paths.betanames = sprintf('%s/BetaNames.txt',paths.results);
paths.sessnums = sprintf('%s/SessNums.txt',paths.results);
paths.combined_cons = sprintf('%s/Betas.nii',paths.results);
paths.thres_file = sprintf('%s/thresholds.mat',paths.results);

paths.realign_params = fullfile(paths.realign,'realign_params.mat');

paths.GLM_noisereg_file = fullfile(paths.GLMdenoise,'GLM_noisereg.mat');
paths.GLM_denoise_params = fullfile(paths.GLMdenoise,'GLM_denoise_params.mat');






%% Stop here if there is no functionnal data
if isempty(SR)
    return
end

%% Check some params
if ~ismember(GLM_params.motion_reg_type,[0 1 2 3])
    error('GLM_params.motion_reg_type value should be 0, 1, 2 or 3')
end

if (n_total_runs < 2) && GLM_params.GLMdenoise && ~paths.in_jack
    error('Cannot run GLMdenoise if the total number of runs is less than n=%i\nChange your parameters to correct this error (GLM_params.GLMdenoise_nPCA_PCs).',2)
end

%% scan info
if exist('scan_info','var')
    clear scan_info scan_log
end
run_id = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        run_id = run_id + 1;
        json_file = [fileparts2(SR(s).filename{r}) '.json'];
        try % fetch all json files and put them into scan_info
            ScanInfos = spm_jsonread(json_file);
            fn = fieldnames(ScanInfos);
            for f = 1:numel(fn)
                field = fn{f};
                scan_info(s,r).(field) = ScanInfos.(field);
            end
        catch ME
            error('.json file of session %02.0f, run %02.0f has a different structure compared to the previous one(s)\nMatlab error:\n%s',SR(s).session,runs(r),getReport(ME,'extended'))
        end
        
        ScanLog_file = fullfile(paths.dataset,'sourcedata',['sub-' paths.subject],sprintf('ses-%02.0f',session),'func',[SR(s).namebase{r} '_ScanLog.mat']);
        
        %         pathstr = fileparts(fileparts(SR(s).filename{r}));
        %         ScanLog_file = fullfile(pathstr,'ScanLogs',[SR(s).namebase{r} '_ScanLog.mat']);
        
        if exist('ScanLog','var'); clear ScanLog; end
        load(ScanLog_file)
        try % fetch all ScanLogs & check if they are all similar
            %             scan_log(s,r) = ScanLog; % this line doesn't work if the
            %             order the fields in the structs is different. Solved with the
            %             four following lines
            fn = fieldnames(ScanLog);
            for f = 1:numel(fn)
                field = fn{f};
                scan_log(s,r).(field) = ScanLog.(field);
            end
            
            if run_id == 1
                ScanLog_ref = ScanLog; % first run is the reference for validation
                nlevels = numel(fieldnames(ScanLog_ref.fMRISTAT));
            else
                if ~strcmp(ScanLog_ref.subject,scan_log(s,r).subject)
                    error('ScanLog of session %02.0f, run %02.0f is for a different subject than the previous one(s)',SR(s).session,runs(r))
                end
                if ~strcmp(ScanLog_ref.TaskName,scan_log(s,r).TaskName)
                    error('ScanLog of session %02.0f, run %02.0f is for a different task than the previous one(s)',SR(s).session,runs(r))
                end
                % try
                %     for L = 1:nlevels
                %         if any(~cellfun(@strcmp,ScanLog_ref.fMRISTAT.(['L' num2str(L)]).names,scan_log(s,r).fMRISTAT.(['L' num2str(L)]).names))
                %             error('Conditions of session %02.0f, run %02.0f, level %i are different from the previous one(s)',SR(s).session,runs(r),L)
                %         end
                %     end
                % catch MEfor
                %     error('Conditions do not match between runs\nMatlab error:\n%s',getReport(MEfor,'extended'))
                % end
            end
        catch ME
            error('ScanLog of session %02.0f, run %02.0f has a different structure compared to the previous one(s)\nMatlab error:\n%s',SR(s).session,runs(r),getReport(ME,'extended'))
        end
    end
end

%%%%%%% Until the dcm2nii converter at CERIMED gives us those info in the json files, read the header using FSL:
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        fname = spm_BIDS(BIDS,'data','sub',paths.subject,'ses',sprintf('%02.0f',session),'run',sprintf('%02.0f',run),'type','bold','task',paths.task);
        [~,comres] = system([paths.FSL_prefix 'fslval ' fname{1} ' dim4']);
        scan_info(s,r).NumberOfVolumesInFile = str2double(comres);
        [~,comres] = system([paths.FSL_prefix 'fslval ' fname{1} ' pixdim1']);
        scan_info(s,r).VoxelSizeX = str2double(comres);
        [~,comres] = system([paths.FSL_prefix 'fslval ' fname{1} ' pixdim2']);
        scan_info(s,r).VoxelSizeY = str2double(comres);
        [~,comres] = system([paths.FSL_prefix 'fslval ' fname{1} ' dim3']);
        scan_info(s,r).NumberOfSlices = str2double(comres);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(unique([scan_info.RepetitionTime])) ~= 1)% || (length(unique([scan_info.nslices])) ~= 1)
    warning('Runs have different TRs')
end
if (length(unique([scan_info.VoxelSizeX])) ~= 1) || (length(unique([scan_info.VoxelSizeY])) ~= 1) || (length(unique([scan_info.SliceThickness])) ~= 1)
    warning('Runs have different voxel sizes')
end

PEDs = cell(n_total_runs,1);
sr = cell(n_total_runs,2);
i = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        i = i + 1;
        PEDs{i} = scan_info(s,r).PhaseEncodingDirection;
        sr{i,1} = session;
        sr{i,2} = run;
    end
end

if n_total_runs > 1
    if ~isequal(PEDs{:})
        warning('Runs have different Phase Encoding Directions:')
        eval('[{''ses'' ''run'' ''PED''};sr PEDs]')
    end
end

%% Display categories at each level
fprintf('\n')
for L = 1:nlevels
    if L == nlevels && L > 1
        fprintf('### Categories in level %i:\n\tAll sounds (sound level)\n',L)
    else
        fprintf('### Categories in level %i:\n\t%s\n',L,strjoin(ScanLog_ref.fMRISTAT.(['L' num2str(L)]).names,', '))
    end
end
fprintf('\n')

%% Contrasts
contrasts.nconds = numel(ScanLog_ref.fMRISTAT.(['L' num2str(contrasts.L)]).names);
contrasts.names = cell(contrasts.nconds,1);
contrasts.names{1} = 'sound_vs_silence';
contrasts.weights = zeros(contrasts.nconds,contrasts.nconds);
contrasts.weights(1,:) = [ones(1,contrasts.nconds-1) -(contrasts.nconds-1)]; % silence always at the end
for c = 2:contrasts.nconds % each condition vs silence
    contrasts.names{c} = [ScanLog_ref.fMRISTAT.(['L' num2str(contrasts.L)]).names{c-1} '_vs_' ScanLog_ref.fMRISTAT.(['L' num2str(contrasts.L)]).names{end}];
    contrasts.weights(c,c-1) = 1;
    contrasts.weights(c,end) = -1;
end

% c = 0; % each condition
% for i = contrasts.nconds+1:contrasts.nconds+contrasts.nconds
%     c = c + 1;
%     contrasts.names{i} = ScanLog_ref.fMRISTAT.(['L' num2str(contrasts.L)]).names{c};
%     contrasts.weights(i,c) = 1;
% end

if exist('contrast_names','var')
    contrast_names = strrep(contrast_names,' ','_'); % deblank
    contrasts.names = [contrasts.names;contrast_names];
end
if exist('contrast_weights','var')
    contrasts.weights = [contrasts.weights;contrast_weights];
end

% Restrict & exclude
if isfield(contrasts,'restrict')
    restrict_L = [contrasts.restrict(:).L];
    restrict_i = find(ismember(restrict_L,contrasts.L));
    if numel(restrict_i) > 1
        error('You can only restrict once at the same level!')
    end
else
    restrict_i = [];
end

if isfield(contrasts,'exclude')
    exclude_L = [contrasts.exclude(:).L];
    exclude_i = find(ismember(exclude_L,contrasts.L));
    if numel(exclude_i) > 1
        error('You can only exclude once at the same level!')
    end
else
    exclude_i = [];
end

if ~isempty(restrict_i) && ~isempty(exclude_i)
    error('You cannot have restrictions AND exclusions at the same level!!')
elseif ~isempty(restrict_i)
    cc = ~contrasts.restrict(restrict_i).c;
elseif ~isempty(exclude_i)
    cc = logical(contrasts.exclude(exclude_i).c);
else
    cc = [];
end


if ~isempty(cc)
    contrasts.nconds = contrasts.nconds - sum(cc);
    for c = 1:size(contrasts.weights,1)
        sc = sum(contrasts.weights(c,:));
        if any(contrasts.weights(c,cc))
            fprintf('Contrast ''%s'' modified by restriction or exclusion\n',contrasts.names{c})
        end
        contrasts.weights(c,cc) = 0;
        neg = sum(contrasts.weights(c,contrasts.weights(c,:) < 0));
        pos = sum(contrasts.weights(c,contrasts.weights(c,:) > 0));
        if ~sc
            if ~neg || ~pos
                contrasts.weights(c,:) = 0;
            elseif neg && pos
                contrasts.weights(c,contrasts.weights(c,:) < 0) = contrasts.weights(c,contrasts.weights(c,:) < 0) .* (pos / -neg);
                contrasts.weights(c,contrasts.weights(c,:)~=0) = contrasts.weights(c,contrasts.weights(c,:)~=0) ./ min(abs(contrasts.weights(c,contrasts.weights(c,:)~=0)));
            end
        end
    end
    contrasts.names(~any(contrasts.weights,2)) = [];
    contrasts.weights(~any(contrasts.weights,2),:) = [];
    contrasts.names
    contrasts.weights
end



if isfield(contrasts,'restrict') || isfield(contrasts,'exclude')
    contrasts.weights(:,end+1) = 0; % add contrast column for exluded or restricted condition
end



%% Realign & Reslice wrap direction
PED = unique([scan_info.PhaseEncodingDirection]);
if strcmp(PED,'i')
    realign_flags.wrap = [1 0 0];
elseif strcmp(PED,'j')
    realign_flags.wrap = [0 1 0];
elseif strcmp(PED,'k')
    realign_flags.wrap = [0 0 1];
end
reslice_flags.wrap = realign_flags.wrap;

%% Use scan info for some parameters
mvt_params.TR = scan_info(1,1).RepetitionTime; % TR (in sec) of the first session first run (all TRs should be the same)
mvt_params.voxel_size(1) = unique([scan_info.VoxelSizeX]); % get voxel size of the protocol
mvt_params.voxel_size(2) = unique([scan_info.VoxelSizeY]);
mvt_params.voxel_size(3) = unique([scan_info.SliceThickness]);
realign_flags.sep = realign_flags.sep * min(mvt_params.voxel_size); %  % from voxels to mm, the default separation (mm) to sample the images.
realign_flags.fwhm = realign_flags.fwhm * min(mvt_params.voxel_size); % % from voxels to mm, The FWHM (mm) applied to the images before estimating the realignment parameters.
mvt_params.max_trans_x = mvt_params.voxel_size(1) * mvt_params.max_trans_vox_x; % from voxels to mm
mvt_params.max_trans_y = mvt_params.voxel_size(2) * mvt_params.max_trans_vox_y;
mvt_params.max_trans_z = mvt_params.voxel_size(3) * mvt_params.max_trans_vox_z;
mvt_params.after_mvt_lag = ceil(mvt_params.after_mvt_lag / (mvt_params.TR)); % from seconds to TRs, number of steady volumes to remove after a movement period

smooth_params.fwhm = smooth_params.fwhm .* mvt_params.voxel_size; % Specify the FWHM of the Gaussian smoothing kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions

if ~isfield(paths,'smoothed_suffix')
    paths.smoothed = sprintf('%s/smoothed_%s',paths.analysis,sprintf('%g',round(smooth_params.fwhm))); % folder where the smoothed resliced motion corrected images will be
else
    paths.smoothed = sprintf('%s/smoothed_%s_%s',paths.analysis,sprintf('%g',round(smooth_params.fwhm)),paths.smoothed_suffix); % folder where the smoothed resliced motion corrected images will be
end


%% write log file parameters
if ~exist(paths.analysis,'dir');mkdir(paths.analysis);end % create folder if non-existant
fid = fopen(paths.log_file_params,'w');
fprintf(fid,'### Subject: %s\n',paths.subject);
fprintf(fid,'### Dataset path: %s\n',paths.dataset);
fprintf(fid,'### Realign name: %s\n',paths.realign_name);

fprintf(fid,'\n\n### Session(s):\n');
fprintf(fid,'%02.0f ',SR.session);

fprintf(fid,'\n### Run(s):\n');
for s = 1:length(SR)
    fprintf(fid,'%02.0f ',SR(s).runs);
    fprintf(fid,'\n');
end

fprintf(fid,'\n### Total number of runs: %i\n',n_total_runs);

fprintf(fid,'\n\n### Realign parameters:\n');
print_struct(realign_flags,fid)
fprintf(fid,'\n\n### Reslice parameters:\n');
print_struct(reslice_flags,fid)
fprintf(fid,'\n\n### Movement parameters:\n');
print_struct(mvt_params,fid)
fprintf(fid,'\n\n### Smoothing parameters:\n');
fprintf(fid,'fwhm in voxels:'); fprintf(fid,' %g',smooth_params.fwhm ./ mvt_params.voxel_size); fprintf(fid,'\n');
print_struct(smooth_params,fid)
fprintf(fid,'\n\n### GLM parameters:\n');
print_struct(GLM_params,fid)
fprintf(fid,'\n\n### Coregistration:\n');
print_struct(coreg_params,fid)
fprintf(fid,'\n\n### Contrasts:\n');
print_struct(contrasts,fid)
fprintf(fid,'\n\n\n### Scan infos of first scan:\n');
print_struct(scan_info(1,1),fid)

fclose(fid);







