fprintf('Preparing to run GLMdenoise...\n')

design = cell(1,n_total_runs);
data = cell(1,n_total_runs);
pca_tissues = {'white';'csf'};

if GLM_params.GLMdenoise_nPCs.(pca_tissues{1}) || GLM_params.GLMdenoise_nPCs.(pca_tissues{2})
    opt.extraregressors = cell(1,n_total_runs);
end
mean_dur = 0;
mean_TR = 0;

volume_correspondance = cell(1,n_total_runs);

opt = [];

run_id = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        run_id = run_id + 1;

        if isfield(scan_log,'Scan_trigs')
            frametimes = scan_log(s,r).Scan_trigs;
        else
            frametimes = 0:scan_info(s,r).RepetitionTime:scan_info(s,r).RepetitionTime * (scan_info(s,r).NumberOfVolumesInFile - 1);
        end
        
        if (max(unique(diff(frametimes))) / scan_info(s,r).RepetitionTime) > 1.1
            variable_TR = 1;
        else
            variable_TR = 0;
        end

        [~,name,ext] = fileparts(SR(s).filename{r});
        smoothed_file = fullfile(paths.smoothed,[smooth_params.prefix reslice_flags.prefix name ext]);

        real_TR = scan_info(s,r).RepetitionTime * 1000; % in msec (which is the precision of the scanner)
        real_TR_sec = real_TR / 1000; % back to sec
        frametimes = reshape(frametimes,1,max(size(frametimes))); % row vector
        onsets = scan_log(s,r).SPM.L1.onsets; % sub levels don't matter here to derive GLMdenoise regressors
        
        if variable_TR
            T1norm_file = fullfile(paths.T1_normalized,['T1norm_' smooth_params.prefix reslice_flags.prefix name ext]);

            jwd = 1.2; % factor by which the TR is multiplied to differentiate silent periods & continuous acquisition
            df = diff(frametimes) * 1000; % sec to msec
            fake_TR = round(mean(df(df<real_TR*jwd))); % round: reduce TR precision to msec
                
            j = 1;
            volume_correspondance{run_id} = zeros(size(frametimes));
            for i = 1:length(frametimes)
                volume_correspondance{run_id}(j) = i;
                if i < length(frametimes)
                    if df(i) < (real_TR * jwd)
                        j = j + 1;
                    else
                        n_fake_vols = ceil((df(i) - fake_TR) / fake_TR);
                        j = j + n_fake_vols + 1;
                    end
                end
            end

            new_frametimes = 0:fake_TR:fake_TR * (length(volume_correspondance{run_id}) - 1);
            new_frametimes = new_frametimes ./ 1000; % back to sec
            fake_TR = fake_TR / 1000; % back to sec
                
            if ~exist(T1norm_file,'file')
                fprintf('Loading session %i, run %i...\n',session,run)
                P = spm_vol(smoothed_file);
                Y = spm_read_vols(P);

                for i = 1:length(scan_log(s,r).Stim_onsets)
                    onset = scan_log(s,r).Stim_onsets(i);
                    corresp_vol = find((frametimes - onset) > 0,1); % first frametime that occured after the stim onset
                    onset_to_frametime = frametimes(corresp_vol) - onset;
                    new_onset = new_frametimes(volume_correspondance{run_id} == corresp_vol) - onset_to_frametime;
                    onsets{scan_log(s,r).SPM.L1.sequence(i)}(onsets{scan_log(s,r).SPM.L1.sequence(i)} == onset) = new_onset; % change to new structure
        %             Stim_onsets(i) = new_onset;
                end

                % load PCA volumes
                [PCA_mat_file,PCA_nii_file] = select_PCA_file(smoothed_file,paths,'no_mask');
                if ~exist(PCA_mat_file,'file')
                    run_PCA_image(smoothed_file,paths,'no_mask',paths.no_mask);
                end

                fmristat_PCA = load(PCA_mat_file); % load fMRISTAT PCA components
                
                Ppca = spm_vol(PCA_nii_file);
                Ypca = spm_read_vols(Ppca);

                PC = 1; % only consider the first PC
                Ypca1 = squeeze(Ypca(:,:,:,PC));

                % find the RMS ratio between min and max RMS for each volume of our time series Y
                maxRMS = 0;
                minRMS = inf;
                for t = 1:size(Y,4)
                    yt = squeeze(Y(:,:,:,t));
                    RMS = sqrt(mean(yt(Ypca1 > 0.99).^2));
                    if RMS > maxRMS
                        maxRMS = RMS;
                    end
                    if RMS < minRMS
                        minRMS = RMS;
                    end
                end
                RMSratio = minRMS / maxRMS;

                % Normalize component with RMSratio (the output w is the factor to apply to the time series in order to temporally normalize it)
                v = -fmristat_PCA.PCs(:,PC);
                v = v + (1 - max(v));
                minv = min(v);
                w = interp1([minv 1],[RMSratio 1],v);

                Ypca(isnan(Ypca)) = 0; % convert NaNs to zeros

                % Apply temporal + spatial normalization (using PCA image)
                for t = 1:size(Y,4)
                    Ypcaw = Ypca(:,:,:,PC);
                    Ypcaw = interp1([1 0 -1],[w(t) 1 (1/w(t))],Ypcaw); % weight PCA image with w
                    % Y(:,:,:,t) = Y(:,:,:,t) .* Ypcaw; % normalize (both spatial & temporal)
                    Y(:,:,:,t) = Y(:,:,:,t) .* w(t); % no spatial normalization (temporal only)
                end
                clear Ypcaw Ypca



    %             fprintf('Writing series...\n') %%%%%%%% MIGHT BE UNNECESSARY TO WRITE THESE VOLUMES
    %             new_P = P;
    %             new_fname = fullfile(paths.T1_normalized,['T1norm_' smooth_params.prefix reslice_flags.prefix name ext]);
    %             if ~exist(paths.T1_normalized,'dir');mkdir(paths.T1_normalized);end % create folder if non-existant
    %             for i = 1:size(Y,4)
    %                 new_P(i) = P(1);
    %                 new_P(i).fname = new_fname;
    %                 new_P(i).n = [i 1];
    %                 spm_write_vol(new_P(i),squeeze(Y(:,:,:,i)));
    %             end


    %             % Kendrick's way of normalizing
    %             mean_before = mean(Y(:));
    %             X = Y;
    %             for t = 1:size(Y,4)
    %                 X(:,:,:,t) = unitlength(Y(:,:,:,t));
    %             end
    %             mean_after = mean(X(:));
    %             X = X * (mean_before / mean_after);


                % create new time series with fake volumes to fill the blanks during silences
                new_Y = zeros(size(Y,1),size(Y,2),size(Y,3),length(new_frametimes));

                new_Y(:,:,:,volume_correspondance{run_id} > 0) = Y; % place "real" volumes
                for i = 1:length(new_frametimes) % add "fake" volumes
                    if volume_correspondance{run_id}(i) == 0
                        I = find(volume_correspondance{run_id}(1:i) > 0,1,'last'); % fake volumes are repetitions of the last volume before the silent gap
                        new_Y(:,:,:,i) = Y(:,:,:,volume_correspondance{run_id}(I));
                    end
                end
                Y = new_Y;
                clear new_Y

                fprintf('Writing T1-normalized series...\n')
                new_P = P;
                if ~exist(paths.T1_normalized,'dir');mkdir(paths.T1_normalized);end % create folder if non-existant
                for i = 1:length(new_frametimes)
                    new_P(i) = P(1);
                    new_P(i).fname = T1norm_file;
                    new_P(i).n = [i 1];
                    spm_write_vol(new_P(i),squeeze(Y(:,:,:,i)));
                end
            else % if ~exist(T1norm_file,'file')
                fprintf('Loading T1norm_file session %i, run %i...\n',session,run)
                P = spm_vol(T1norm_file);
                Y = spm_read_vols(P);
            end % if ~exist(T1norm_file,'file')
            
            % include PCA PCs
            j = 0;
            clear confounds
            for i = 1:numel(pca_tissues)
                if GLM_params.GLMdenoise_nPCs.(pca_tissues{i}) > 0
                    PCA_mat_file = select_PCA_file(T1norm_file,paths,pca_tissues{i});
                    if ~exist(PCA_mat_file,'file')
                        run_PCA_image(T1norm_file,paths,pca_tissues{i},paths.func.pca_mask.(pca_tissues{i}));
                    end
                    fmristat_PCA.(pca_tissues{i}) = load(PCA_mat_file);
                    
                    if GLM_params.GLMdenoise_nPCs.(pca_tissues{i}) <= size(fmristat_PCA.(pca_tissues{i}).PCs,2)
                        nPCs.(pca_tissues{i}) = GLM_params.GLMdenoise_nPCs.(pca_tissues{i});
                    else
                        error('The wanted number of PCA PCs (GLM_params.nPCs.%s) is greater than the number of PCs available (n = %i)\nChange your parameters to correct this error.',pca_tissues{i},size(fmristat_PCA.PCs,2))
                    end
                    if ~exist('confounds','var')
                        confounds = fmristat_PCA.(pca_tissues{i}).PCs(:,1:nPCs.(pca_tissues{i}));
                    else
                        confounds(:,j+1:j+nPCs.(pca_tissues{i})) = fmristat_PCA.(pca_tissues{i}).PCs(:,1:nPCs.(pca_tissues{i}));
                    end
                    j = j + nPCs.(pca_tissues{i});
                end
            end
            if exist('confounds','var')
                opt.extraregressors{run_id} = confounds; % add fMRISTAT PCA components as regressors
            end
            
            
            mean_TR = mean_TR + fake_TR;
        else % if not variable_TR
            fprintf('Loading session %i, run %i...\n',session,run)
            P = spm_vol(smoothed_file);
            Y = spm_read_vols(P);
            
            % include PCA PCs
            j = 0;
            clear confounds
            for i = 1:numel(pca_tissues)
                if GLM_params.GLMdenoise_nPCs.(pca_tissues{i}) > 0
                    PCA_mat_file = select_PCA_file(smoothed_file,paths,pca_tissues{i});
                    if ~exist(PCA_mat_file,'file')
                        run_PCA_image(smoothed_file,paths,pca_tissues{i},paths.func.pca_mask.(pca_tissues{i}));
                    end
                    fmristat_PCA.(pca_tissues{i}) = load(PCA_mat_file);
                    
                    if GLM_params.GLMdenoise_nPCs.(pca_tissues{i}) <= size(fmristat_PCA.(pca_tissues{i}).PCs,2)
                        nPCs.(pca_tissues{i}) = GLM_params.GLMdenoise_nPCs.(pca_tissues{i});
                    else
                        error('The wanted number of PCA PCs (GLM_params.nPCs.%s) is greater than the number of PCs available (n = %i)\nChange your parameters to correct this error.',pca_tissues{i},size(fmristat_PCA.PCs,2))
                    end
                    if ~exist('confounds','var')
                        confounds = fmristat_PCA.(pca_tissues{i}).PCs(:,1:nPCs.(pca_tissues{i}));
                    else
                        confounds(:,j+1:j+nPCs.(pca_tissues{i})) = fmristat_PCA.(pca_tissues{i}).PCs(:,1:nPCs.(pca_tissues{i}));
                    end
                    j = j + nPCs.(pca_tissues{i});
                end
            end
            if exist('confounds','var')
                opt.extraregressors{run_id} = confounds; % add fMRISTAT PCA components as regressors
            end


            mean_TR = mean_TR + real_TR_sec;
        end % if variable_TR
    	  
        % design{run_id} = onsets';
        design{run_id} = onsets(1:end-1)'; % don't include silent condition
        data{run_id} = single(Y);

        mean_dur = mean_dur + mean(scan_log(s,r).Stim_durations);
    end
end

mean_dur = mean_dur / n_total_runs;
mean_TR = mean_TR / n_total_runs;


% opt.numpcstotry = 20;
% opt.numboots = 0;

if ~exist(paths.results,'dir');mkdir(paths.results);end % create folder if non-existant
if ~exist(paths.GLMdenoise,'dir');mkdir(paths.GLMdenoise);end % create folder if non-existant

% hrf_canon = getcanonicalhrf(mean_dur,mean_TR)';
% time=(0:length(hrf_canon)-1);
% hrf_struct=fmridesign(time,0,[1 0],[],GLM_params.hrf);
% hrf = squeeze(hrf_struct.X(:,1,1,1));
% if strcmpi(paths.contrast_agent,'MION')
%     hrf = -hrf;
% end

% GLM_denoise_results = GLMdenoisedata(design,data,mean_dur,mean_TR,'assume',hrf,opt,fullfile(paths.GLMdenoise,'figures')); % [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)

GLM_denoise_results = GLMdenoisedata(design,data,mean_dur,mean_TR,'assume',[],opt,fullfile(paths.GLMdenoise,'figures')); % [results,denoiseddata] = GLMdenoisedata(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,figuredir)

GLM_noisereg = cellfun(@(x) x(:,1:GLM_denoise_results.pcnum),GLM_denoise_results.pcregressors,'UniformOutput',0);


if variable_TR
    for i = 1:n_total_runs
        GLM_noisereg{i} = GLM_noisereg{i}(volume_correspondance{i}>0,:);
    end
end

save(paths.GLM_noisereg_file,'GLM_noisereg')

% filename = fullfile(paths.GLMdenoise,'GLM_denoise_results.mat');
% save(filename,'GLM_denoise_results')

save(paths.GLM_denoise_params,'SR')

% plot(Stim_onsets,[ones(size(Stim_onsets))],'o')
% hold on
% plot(new_frametimes,volume_correspondance>0)

