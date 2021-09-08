%% Run GLMdenoise if needed
if GLM_params.GLMdenoise
    runit = 1;
    if exist(paths.GLM_denoise_params,'file')
        tocheck = load(paths.GLM_denoise_params);
        % if isequal(SR,tocheck.SR) % && isequal(sessions,tocheck.sessions) % || paths.in_jack
        %     runit = 0;
        % end
        if all(ismember([SR.session],[tocheck.SR.session]))
            runit = 0;
            for i = 1:numel(tocheck.SR)
                srs = find([SR.session] == tocheck.SR(i).session);
                if ~isempty(srs)
                    if ~all(ismember(SR(srs).runs,tocheck.SR(i).runs))
                        runit = 1;
                    end
                end
            end
        end
    end
    
    if runit
        run_GLMdenoise
    end
end

%% Global GLM parameters
% ncons = numel(contrasts.names);

if GLM_params.AR > 1
    which_stats = '_mag_t _mag_ef _mag_sd _cor _AR';
else
    which_stats = '_mag_t _mag_ef _mag_sd';
end

fwhm_cor = Inf; % assumes the frames are uncorrelated
n_trends = [GLM_params.n_temporal_trends 1 1];

%% Check if this particular analysis (GLM parameters) already exists
paths.results_indiv_runs = fullfile(paths.results,'individual_trials'); % folder where the results for individual runs will be
paths.GLM_params = fullfile(paths.results_indiv_runs,'GLM_params.mat');

same_analysis = 0;
if exist(paths.GLM_params,'file') == 2
    load(paths.GLM_params)
    if isequaln(GLM_params,saved_GLM_params)
        same_analysis = 1;
    end
end

%% Files
if ~exist(paths.results,'dir');mkdir(paths.results);end % create folder if non-existant
if ~exist(paths.results_indiv_runs,'dir');mkdir(paths.results_indiv_runs);end % create folder if non-existant
out_name = sprintf('sub-%s_task-%s_res-%s',paths.subject,paths.task,paths.results_name);
labels_file = fullfile(paths.results_indiv_runs,[out_name '_Labels.txt']);

if ~same_analysis

    %% Labels & categories txt files
    labels_fid = fopen(labels_file,'w');
    categ_fid = fopen(fullfile(paths.results_indiv_runs,[out_name '_Categories.txt']),'w');
    fprintf(labels_fid,'trial\t');
    Lnames = fieldnames(scan_log(1,1).fMRISTAT);
    nrows = 0;
    for i = 1:length(Lnames)
        fprintf(labels_fid,'%s_num\t',Lnames{i});
        fprintf(labels_fid,'%s\t',Lnames{i});
        fprintf(categ_fid,Lnames{i});
        if i < length(Lnames)
            fprintf(categ_fid,'\t');
        else
            fprintf(categ_fid,'\n');
        end
        if numel(scan_log(1,1).fMRISTAT.(['L' num2str(i)]).names) > nrows
            nrows = numel(scan_log(1,1).fMRISTAT.(['L' num2str(i)]).names);
        end
    end
    fprintf(labels_fid,'session\trun\trunID\n');


    for ro = 1:nrows-1 % skip silent condition
        for i = 1:length(Lnames)
            if ro < numel(scan_log(1,1).fMRISTAT.(['L' num2str(i)]).names)
                fprintf(categ_fid,scan_log(1,1).fMRISTAT.(['L' num2str(i)]).names{ro});
            end
            if i < length(Lnames)
                fprintf(categ_fid,'\t');
            end
        end
        if ro < nrows
            fprintf(categ_fid,'\n');
        end
    end
    fclose(categ_fid);




    %% Analyze each run
    trial_num = 0;
    run_id = 0;
    silent_cond = find(strcmp(scan_log(1,1).fMRISTAT.(['L' num2str(contrasts.L)]).names,'silence'));
    for s = 1:numel(SR)
        session = SR(s).session;
        runs = SR(s).runs;
        for r = 1:length(runs)
            run = runs(r);

            run_id = run_id + 1;

            [~,name,ext] = fileparts(SR(s).filename{r});
            smoothed_file = fullfile(paths.smoothed,[smooth_params.prefix reslice_flags.prefix name ext]);

            if isfield(scan_log(s,r),'Scan_trigs') && ~isempty(scan_log(s,r).Scan_trigs)
                frametimes = scan_log(s,r).Scan_trigs;
            else
                frametimes = 0:scan_info(s,r).RepetitionTime:scan_info(s,r).RepetitionTime * (scan_info(s,r).NumberOfVolumesInFile - 1);
            end
            frametimes = reshape(frametimes,1,max(size(frametimes))); % row vector
            slicetimes = scan_info(s,r).SliceTiming;
            slicetimes = reshape(slicetimes,1,max(size(slicetimes))); % row vector
            events = scan_log(s,r).fMRISTAT.(['L' num2str(contrasts.L)]).events; % the design matrix

            %% Write labels file & make contrasts
            contrasts_names = cell(length(events(events(:,1) ~= silent_cond)),1);
            this_run_i = 0;
            for i = 1:size(events,1)
                not_silence = 0;
                for L = 1:length(Lnames)
                    etype = scan_log(s,r).fMRISTAT.(Lnames{L}).events(i,1);
                    ename = scan_log(s,r).fMRISTAT.(Lnames{L}).names{etype};
                    if not_silence && L > 1 && strcmp(ename,'silence')
                        keyboard
                    end
                    if ~strcmp(ename,'silence')
                        if L == 1
                            not_silence = 1;
                            trial_num = trial_num + 1;
                            this_run_i = this_run_i + 1;
                            fprintf(labels_fid,'%i\t',trial_num);
                        end
                        fprintf(labels_fid,'%i\t%s\t',etype,ename);
                    end
                end
                if not_silence
                    fprintf(labels_fid,'%i\t%i\t%i\n',session,run,run_id);
                    contrasts_names{this_run_i,1} = sprintf('trial-%05.0f',trial_num);
                    events(i,1) = this_run_i; % number of this trial which is not a silence
                else
                    events(i,1) = inf; % put inf when silence
                end
            end

            % change contrasts weights to match this design
            contrasts_weights = eye(this_run_i);
            contrasts_weights(:,end+1) = ones(this_run_i,1) .* -1;

            % change events to set silence condition
            events(events(:,1)==inf,1) = this_run_i + 1;

            %% Perform GLM
            ncons_this_run = numel(contrasts_names);

            if ncons_this_run > 0

                out_base = cell(ncons_this_run,1);

                for con = 1:ncons_this_run
                    out_base{con} = fullfile(paths.results_indiv_runs,[out_name '_' contrasts_names{con}]);
                end

                fprintf('\nPerforming GLM on session %i, run %i:\n',session,run)
                out_base = char(out_base);

                X_cache = fmridesign(frametimes,slicetimes,events,[],GLM_params.hrf);

                
                % if PCA PCs have to be passed as regressors
                for i = 1:numel(pca_tissues)
                    if GLM_params.nPCs.(pca_tissues{i}) > 0
                        PCA_mat_file = select_PCA_file(smoothed_file,paths,pca_tissues{i});
                        if ~exist(PCA_mat_file,'file')
                            run_PCA_image(smoothed_file,paths,pca_tissues{i},paths.func.pca_mask.(pca_tissues{i}));
                        end
                        fmristat_PCA.(pca_tissues{i}) = load(PCA_mat_file);
                        
                        if GLM_params.nPCs.(pca_tissues{i}) <= size(fmristat_PCA.(pca_tissues{i}).PCs,2)
                            nPCs.(pca_tissues{i}) = GLM_params.nPCs.(pca_tissues{i});
                        else
                            error('The wanted number of PCA PCs (GLM_params.nPCs.%s) is greater than the number of PCs available (n = %i)\nChange your parameters to correct this error.',pca_tissues{i},size(fmristat_PCA.PCs,2))
                        end
                    else
                        nPCs.(pca_tissues{i}) = 0;
                    end
                end



                if GLM_params.motion_reg_type == 0
                    nMR = 0;
                elseif GLM_params.motion_reg_type == 1
                    nMR = 6;
                elseif GLM_params.motion_reg_type == 2
                    nMR = 12;
                elseif GLM_params.motion_reg_type == 3
                    nMR = 24;
                end
                if GLM_params.GLMdenoise % if using GLM denoise to derive noise regressors
                    load(paths.GLM_noisereg_file);
                    nGLMdenoisePCs = size(GLM_noisereg{run_id},2);
                else
                    nGLMdenoisePCs = 0;
                end


                RealignRef_file = fullfile(paths.realign,[name '_RealignRef.mat']);
                load(RealignRef_file)
                if GLM_params.pass_rejected_vols_as_regressors
                    confounds = zeros(size(selected_vols,1),sum(~selected_vols) + nMR + nPCs.white + nPCs.csf + nGLMdenoisePCs); % all regressors
                else
                    confounds = zeros(size(selected_vols,1),nMR + nPCs.white + nPCs.csf + nGLMdenoisePCs); % all regressors
                end

                if ~isempty(confounds)
                    j = 0;
                    if GLM_params.pass_rejected_vols_as_regressors
                        if GLM_params.all_trials_are_ok
                            % consider all volumes within a trial as ok
                            for i = 1:length(scan_log(s,r).Stim_onsets)
                                selected_vols((scan_log(s,r).Scan_trigs > scan_log(s,r).Stim_onsets(i)) & (scan_log(s,r).Scan_trigs < scan_log(s,r).Reward_onsets(i))) = 1;
                            end
                        end
                        for i = 1:length(selected_vols)
                            if ~selected_vols(i)
                                j = j + 1;
                                confounds(i,j) = 1;
                            end
                        end
                    end

                    % motion regressors
                    if nMR > 0
                        rp_file = fullfile(paths.rp_files,['rp_' name '.txt']);
                        Q = load(rp_file);
                        q6 = j+1:j+6;
                        confounds(:,q6) = Q;
                        j = j + 6;
                        if nMR > 6
                            confounds(:,j+1:j+6) = confounds(:,q6) .^ 2; % squared
                            j = j + 6;
                            if nMR > 12
                                confounds(2:end,j+1:j+6) = diff(confounds(:,q6)); % derivate
                                confounds(2:end,j+7:j+12) = diff(confounds(:,q6) .^ 2); % derivate of squared
                                j = j + 12;
                            end
                        end
                    end


                    % white matter & csf PCA regressors
                    for i = 1:numel(pca_tissues)
                        if nPCs.(pca_tissues{i}) > 0
                            confounds(:,j+1:j+nPCs.(pca_tissues{i})) = fmristat_PCA.(pca_tissues{i}).PCs(:,1:nPCs.(pca_tissues{i}));
                            j = j + nPCs.(pca_tissues{i});
                        end
                    end

                    % GLMdenoise regressors
                    if GLM_params.GLMdenoise
                        confounds(:,j+1:j+nGLMdenoisePCs) = GLM_noisereg{run_id};
                    end
                end

                if ~exist(paths.results,'dir');mkdir(paths.results);end % create folder if non-existant
                if ~exist(paths.results_indiv_runs,'dir');mkdir(paths.results_indiv_runs);end % create folder if non-existant

                % run GLM
                df = fmrilm(smoothed_file,out_base,X_cache,contrasts_weights,[],which_stats,fwhm_cor,n_trends,confounds,[],[],[],GLM_params.AR);
                if isempty(df)
                    error('Non-estimable contrast in FMRISTAT, session %i, run %i might be due to too many confounds',session,run)
                end
            end
        end
    end

    fclose(labels_fid);

    %% Save GLM & contrast parameters
    saved_GLM_params = GLM_params;
    if ~paths.in_jack
        save(paths.GLM_params,'saved_GLM_params')
    end
end % if ~same_analysis

%% Write log file
if ~exist(paths.results_multi,'dir');mkdir(paths.results_multi);end % create folder if non-existant
log_file_results = sprintf('%s/Results_%s.log',paths.results_multi,paths.results_name);
fid = fopen(log_file_results,'w');
fprintf(fid,'Total number of runs analyzed: %i\n\n',n_total_runs);
% fprintf(fid,'Number of events per condition:\n');
% a = char(scan_log(1,1).fMRISTAT.(['L' num2str(contrasts.L)]).names);
% for c = 1:size(a,1)
%     fprintf(fid,'%s %i\n',a(c,:),n_events_per_cond(c));
% end

%% Get conditions for each trial
labels = tdfread(labels_file);
this_L_conds = labels.(['L' num2str(contrasts.L) '_num']);

%% Input files
dir_ef = dir(fullfile(paths.results_indiv_runs,'*_mag_ef.nii'));
dir_sd = dir(fullfile(paths.results_indiv_runs,'*_mag_sd.nii'));
ntrial = length(dir_ef);
if ntrial ~= length(this_L_conds)
    error('Numbers of input files & number of trials mismatch!')
end
infiles_ef = cell(ntrial,1);
infiles_sd = cell(ntrial,1);
for i = 1:ntrial
    infiles_ef{i} = fullfile(paths.results_indiv_runs,dir_ef(i).name);
    infiles_sd{i} = fullfile(paths.results_indiv_runs,dir_sd(i).name);
end
infiles_ef = char(infiles_ef);
infiles_sd = char(infiles_sd);
        
%% Combine trials fixed effects for each contrast %%%
MS_which_stats = '_t _ef _sd';
ws = textscan(MS_which_stats,'%s','Delimiter',' ');

% p_val_peak = GLM_params.p_val_peak;
setenv('ANTSPATH',paths.ANTS_path);
ncons = numel(contrasts.names);
for con = 1:ncons
    if ~(paths.in_jack && ~strcmp(contrasts.names{con},'sound_vs_silence'))

        % Create design matrix
        if isfield(contrasts,'restrict') || isfield(contrasts,'exclude')
            to_remove = 2;
        else
            to_remove = 1;
        end
        MS_contrast = contrasts.weights(con,1:end-to_remove);
        X = zeros(length(this_L_conds),length(MS_contrast));
        indices = sub2ind(size(X),1:length(this_L_conds),this_L_conds');
        X(indices) = 1;
        
        % Combine inter level conditions (restrict to conditions of other levels)
        if isfield(contrasts,'restrict')
            for res = 1:length(contrasts.restrict)
                ok_conds = find(contrasts.restrict(res).c);
                ok_trials = ismember(labels.(['L' num2str(contrasts.restrict(res).L) '_num']),ok_conds);
                X = X .* repmat(ok_trials,1,size(X,2));
            end
        elseif isfield(contrasts,'exclude')
            for res = 1:length(contrasts.exclude)
                ok_conds = find(~contrasts.exclude(res).c);
                ok_trials = ismember(labels.(['L' num2str(contrasts.restrict(res).L) '_num']),ok_conds);
                X = X .* repmat(ok_trials,1,size(X,2));
            end
        end
        
        
        if any((MS_contrast~=0) & any(X)) % if conditions & contrast match
            % Perform multistat
            out_base = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con}]);

            df = multistat(infiles_ef,infiles_sd,[],[],X,MS_contrast,out_base,MS_which_stats,Inf); %tata!
            % Last argument of multistat, FWHM_VARATIO:
            %  - 0 will do no smoothing, and give a purely random effects analysis;
            %  - Inf will do complete smoothing to a global ratio of one, giving a purely fixed effects analysis.

            % find threshold
            [volume, number_voxels] = mask_vol(paths.func_mask); % "As input this function requires a mask_file which is the original Motion Corrected fMRI data"
            [peak_threshold,extent_threshold_vol,Cluster_threshold] = stat_threshold(volume,number_voxels,realign_flags.fwhm,df.t,p_val_peak);
            extent_threshold_vox = round(extent_threshold_vol / abs(prod([scan_info(1,1).VoxelSizeX scan_info(1,1).VoxelSizeY scan_info(1,1).SliceThickness])));

            save([out_base '_threshold.mat'],'peak_threshold','extent_threshold_vol','extent_threshold_vox','Cluster_threshold','df','p_val_peak')


            % Mask output images
            for st = 1:numel(ws{1})
                maths_in_file = sprintf('%s%s.nii',out_base,ws{1}{st});
                if exist([maths_in_file '.gz'],'file') == 2; delete([maths_in_file '.gz']);end
                system(sprintf('%sfslmaths %s -mas %s %s',paths.FSL_prefix,maths_in_file,paths.func_mask,maths_in_file));
                if exist(maths_in_file,'file') == 2; delete(maths_in_file);end
                %                 system(sprintf('gunzip %s',maths_in_file));
            end

            % Get activation
            t_map_file = sprintf('%s_t.nii.gz',out_base);
            ef_map_file = sprintf('%s_ef.nii.gz',out_base);
            sd_map_file = sprintf('%s_sd.nii.gz',out_base);
            [n_vox_mask,n_vox_maskP,meanT,meanTP,maxT] = get_activation(t_map_file,contrasts.names{con},p_val_peak(1),p_val_peak,peak_threshold,extent_threshold_vox,Cluster_threshold,paths);

            fprintf(fid,'\n### AC activation for %s:\n',contrasts.names{con});
            fprintf(fid,'Spatial extent (total) = %i\n',n_vox_mask);
            fprintf(fid,'Spatial extent (peaks) = %i\n',n_vox_maskP);
            fprintf(fid,'Mean t-value (total) = %g\n',meanT);
            fprintf(fid,'Mean t-value (peaks) = %g\n',meanTP);
            fprintf(fid,'Max t-value = %g\n',maxT);

            if ~paths.in_jack
                % Apply transforms to register to anat & template
                fprintf('Applying spatial transforms...\n\n')
                if strcmp(coreg_params.method,'bbr')
                    t_average_to_anat_warp = '';
                elseif strcmp(coreg_params.method,'ants')
                    t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
                else
                    error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
                end
      
                out_file = fullfile(paths.results_multi,['In-' T1w_files.(paths.template_name).name '_sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_t.nii.gz']);
                system(sprintf('%s -i %s -r %s -o %s%s -t %s -n Linear',...
                    fullfile('$ANTSPATH','antsApplyTransforms'),...
                    t_map_file,paths.anat.full,out_file,...
                    t_average_to_anat_warp,paths.average_to_anat_xfm));
             
                % maps to template
                for st = 1:numel(ws{1})
                    in_file = sprintf('%s%s.nii.gz',out_base,ws{1}{st});
                    out_file = fullfile(paths.results_multi,['In-' paths.template_name '_sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} ws{1}{st} '.nii.gz']);
                    system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n Linear',...
                        fullfile('$ANTSPATH','antsApplyTransforms'),...
                        in_file,paths.anat.full,out_file,...
                        paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
                        t_average_to_anat_warp,paths.average_to_anat_xfm));
                end
            end
        end
    end
end

fprintf(fid,'\n### Thresholds:\n');
for i = 1:length(p_val_peak)
    fprintf(fid,'p = %g: Peak threshold = %g\n',p_val_peak(i),peak_threshold(i));
end
fprintf(fid,'Cluster threshold = %g\n',Cluster_threshold);
fprintf(fid,'Extent threshold = %g mm^3\n',extent_threshold_vol);
fprintf(fid,'Extent threshold = %g voxels (in the functional space)\n',extent_threshold_vox);


fclose(fid);



