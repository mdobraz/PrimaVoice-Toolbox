%% Path things
if ~exist(paths.tmp_nii,'dir');mkdir(paths.tmp_nii);end % create tmp folder if non-existant
fid = fopen(paths.log_file_realign,'w');

if ~exist(paths.resliced,'dir');mkdir(paths.resliced);end % create tmp folder if non-existant

%% Set ANTS environnement
setenv('ANTSPATH',paths.ANTS_path);

%% Check if this particular realignment already exists
same_analysis = 0;
if exist(paths.realign_params,'file') == 2
    load(paths.realign_params)
    if isequaln(realign_FSL,saved_realign_FSL) && isequaln(mvt_params,saved_mvt_params)
        same_analysis = 1;
    end
end


%% Find representative scans of each run
ref_scans = cell(n_total_runs,2);
EPIs = cell(n_total_runs,1);
run_count = 0;
vols_selected = zeros(n_total_runs,5);
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    SR(s).ref_volume = nan(length(runs),1);
    SR(s).selected_vols = cell(length(runs),1);
    
    for r = 1:length(runs)
        run = runs(r);
        %% Filename
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);
        if strcmp(ext,'.gz'); [~,bold_name,ext] = fileparts(bold_name); end
        
        %% Check if this particular realignment has already been saved
        RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
        if ~(exist(RealignRef_file,'file') == 2) || ~same_analysis % if the file does not exist, perform realignment estimate
            %% Find ref vol
            fprintf('\nFinding representative volume for session %i, run %i\n',session,run)
            [ref_volume, selected_vols, Q, i1, i2] = find_refvol_FSL(bold_file,realign_FSL,mvt_params,paths);
            
            mcRMS.trans = norm(sqrt(mean(Q(:,1:3).^2))); % norm of RMS of translations
            mcRMS.rot = norm(sqrt(mean(Q(:,4:6).^2))); % norm of RMS of rotations
            
            %% Generate figure
            figureprep([0 0 1280 960]);
            realign_plot_Q(Q,sprintf('Session %i, run %i - Relative to first volume',session,run))
            hold on
            plot([ref_volume ref_volume],[min(min(Q(:,4:6))) max(max(Q(:,4:6)))])
            text(ref_volume+1,(max(max(Q(:,4:6))) / 2),sprintf('representative\nvolume'))
            plot(~selected_vols .* max(max(Q(:,4:6))),'LineWidth',2)
            plot([i1 i2],[(max(max(Q(:,4:6))) * 0.75) (max(max(Q(:,4:6))) * 0.75)])
            hold off
            if ~exist(paths.realign,'dir');mkdir(paths.realign);end % create folder if non-existant
            fig_prefix = fullfile(paths.realign,[bold_name '_Realignment_pre']);
            figurewrite(fig_prefix,[],0,paths.realign); % the 0 is to force eps figure
            
            %% Save realignment files for the run
            save(RealignRef_file,'ref_volume','selected_vols','mcRMS')
        else
            load(RealignRef_file)
        end
        
        %% Write log file
        fprintf(fid,'Session %i, run %i\n',session,run);
        fprintf(fid,'%i / %i volumes selected, reference volume: %i\n',sum(selected_vols),length(selected_vols),ref_volume);
        
        %% Store ref scans
        SR(s).ref_volume(r) = ref_volume;
        SR(s).selected_vols{r} = selected_vols;
        run_count = run_count + 1;
        ref_scans{run_count,1} = bold_file;
        ref_scans{run_count,2} = ref_volume;
        EPIs{run_count} = bold_file;
        
        %% Save vols selected
        if exist('vols_selected','var'); clear vols_selected; end
        vols_selected(run_count,1) = session;
        vols_selected(run_count,2) = run;
        vols_selected(run_count,3) = sum(selected_vols);
        vols_selected(run_count,4) = length(selected_vols) - sum(selected_vols);
        vols_selected(run_count,5) = length(selected_vols);
    end
end
ref_scans_orig = ref_scans;

%% Save vols selected
vols_sel_file = fullfile(paths.realign,['sub-' paths.subject '_task-' paths.task '_VolsSelected']);
save(vols_sel_file,'vols_selected')

saved_realign_FSL = realign_FSL;
saved_mvt_params = mvt_params;
save(paths.realign_params,'saved_realign_FSL','saved_mvt_params')


%% Calculate VDM if fmap present
spm('defaults', 'FMRI');
spm_jobman('initcfg'); % initialization
spm_get_defaults('cmdline',true)
if reslice_flags.unwarp
    has_fmap = zeros(run_count,1);
    all_VDMs = cell(run_count,1);
    run_count = 0;
    for s = 1:numel(SR)
        session = SR(s).session;
        runs = SR(s).runs;
%         if exist('fmap','var'); clear fmap; end
        fieldmaps = spm_BIDS(BIDS,'data','sub',paths.subject,'modality','fmap','type','fieldmap','ses',sprintf('%02.0f',session));
        if ~isempty(fieldmaps)
%             SR(s).fieldmap = cell(length(runs),1);
            if ~exist(paths.fieldmaps,'dir');mkdir(paths.fieldmaps);end % create folder if non-existant
            IntendedFor = cell(numel(fieldmaps),1);
            for fmi = 1:numel(fieldmaps)
                fmap_file = fieldmaps{fmi};
                [fmap_pathstr,fmap_name] = fileparts(fmap_file);
                fmap_json = spm_jsonread(fullfile(fmap_pathstr,[fmap_name '.json']));
                c = ~cellfun(@isempty,strfind(fmap_json.IntendedFor,sprintf('sub-%s_ses-%02.0f_task-%s',paths.subject,session,paths.task)));
                if any(c) % if fieldmap is intended for runs of the current analysis & current session
                    % Intended for
                    IntendedFor{fmi} = fmap_json.IntendedFor(c);
                    % magnitude file (we have to do all this because spm_BIDS does not see the magnitude file)
                    mag_name = [fmap_name(1:strfind(fmap_name,'fieldmap')-1) 'magnitude.nii'];
                    mag_file = fullfile(fmap_pathstr,mag_name);
                    a = dir([mag_file '*']);
                    [~,~,ext] = fileparts(a.name);
                    if strcmp(ext,'.gz'); system(sprintf('gunzip %s',mag_file)); end
                    
                    for epi = 1:numel(IntendedFor{fmi})
                        % echo times (might not be important)
                        [~,epi_name,ext] = fileparts(IntendedFor{fmi}{epi});
                        if strcmp(ext,'.gz'); [~,epi_name] = fileparts(epi_name); end

                        epi_run = str2double(epi_name(strfind(epi_name,'_run-')+5:strfind(epi_name,'_run-')+6));
                        if ismember(epi_run,runs)
                            bold_files = spm_BIDS(BIDS,'data','sub',paths.subject,'ses',sprintf('%02.0f',session),'task',paths.task,'type','bold');
                            epi_pathstr = fileparts(bold_files{1});
                            epi_json_file = fullfile(epi_pathstr,[epi_name '.json']);
                            if exist(epi_json_file,'file') % if run exists
                                epi_json = spm_jsonread(epi_json_file);
                                echo_times = sort([epi_json.EchoTime fmap_json.EchoTime] .* 1000);
                                
                                % Blip direction of EPI & total readout time
                                PED = epi_json.PhaseEncodingDirection;
                                if isempty(strfind(PED,'-'))
                                    blipdir = 1;
                                else
                                    blipdir = -1;
                                end
                                tert = epi_json.TotalReadoutTime * 1000;
                                
                                
                                rv = SR(s).ref_volume(runs == epi_run); % ref volume of this run
                                epi_file = SR(s).filename{SR(s).runs == epi_run};
                                ref_epi = fullfile(paths.tmp_nii,sprintf('%s_%i.nii',epi_name,rv));
                                if exist(ref_epi,'file') == 2; delete(ref_epi);end
                                system(sprintf('%sfslroi %s %s %i 1',paths.FSL_prefix,epi_file,ref_epi,rv-1)); % fslroi indexing starts at 0, hence the '-1' here
                                system(sprintf('gunzip %s',ref_epi));
                                
                                
                                if exist('matlabbatch','var'); clear matlabbatch; end
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.precalcfieldmap = {[fmap_file ',1']};
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.precalcfieldmap.magfieldmap = {[mag_file ',1']};
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = echo_times;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = tert;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 1;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(fileparts(which('spm.m')),'toolbox','FieldMap','T1.nii')};
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = {[ref_epi ',1']};
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'VDMrun';
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
                                matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
                                
                                output_list = spm_jobman('run', matlabbatch); % Run job
                                
                                movefile(output_list{1}.vdmfile{1}{1},paths.fieldmaps);

                                [~,vdmname,vdmext] = fileparts(output_list{1}.vdmfile{1}{1});
                                SR(s).vdmfile{SR(s).runs == epi_run} = fullfile(paths.fieldmaps,[vdmname vdmext]);
                                
                                % Save figures of warped overlaid with VDM & unwarped
                                view_slice_overlay(ref_epi,fullfile(paths.tmp_nii,sprintf('u%s_%i.nii',epi_name,rv)),0,[],1,[],'gray')
                                figurewrite(fullfile(paths.fieldmaps,[epi_name '_unwarped'])) % Using GLMdenoise function
                                view_slice_overlay(ref_epi,SR(s).vdmfile{SR(s).runs == epi_run},0,[],0.2)
                                figurewrite(fullfile(paths.fieldmaps,[epi_name '_VDM'])) % Using GLMdenoise function
                            end
                        end
                    end
                end
            end
            
            for r = 1:length(runs)
                % keep a record of which run has a fieldmap associated with
                run = runs(r);
                run_count = run_count + 1;
                found_run = 0;
                for ci = 1:numel(IntendedFor)
                    c = ~cellfun(@isempty,strfind(IntendedFor{ci},sprintf('sub-%s_ses-%02.0f_task-%s_run-%02.0f',paths.subject,session,paths.task,run)));
                    if any(c)
                        found_run = 1;
                        break
                    end
                end
                if found_run
                    has_fmap(run_count) = 1;
                    all_VDMs{run_count} = SR(s).vdmfile{r};
                end
            end
        end
    end
    
    
    if sum(has_fmap) >= (length(has_fmap) / 2) % there is a majority of runs with a fieldmap
        master_has_fmap = 1;
        VDMs = all_VDMs(logical(has_fmap));
    else % there is a majority of runs without fieldmap
        master_has_fmap = 0;
    end
end


%% Apply VDM if wanted & if fieldmap is present
if reslice_flags.unwarp
    if any(has_fmap)
        if exist('matlabbatch','var'); clear matlabbatch; end
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [reslice_flags.which reslice_flags.mean];
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = reslice_flags.interp;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.wrap = reslice_flags.wrap;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.mask = reslice_flags.mask;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'unwarped_';
        
        run_count = 0;
        for s = 1:numel(SR)
            session = SR(s).session;
            runs = SR(s).runs;
            for r = 1:length(runs)
                run = runs(r);
                run_count = run_count + 1;
                bold_file = SR(s).filename{r};
                [bold_path,bold_name,ext] = fileparts(bold_file);
                out_file = fullfile(bold_path,['unwarped_' bold_name ext]);
                if has_fmap(run_count)
                    Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
                    for vol = 1:numel(Scans)
                        Scans{vol,1} = sprintf('%s,%i',bold_file,vol);
                    end
                    matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.scans = Scans;
                    matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.vdmfile = SR(s).vdmfile(r);
                    % Run job
                    if ~(exist(out_file,'file') == 2) || ~same_analysis
                        spm_jobman('run', matlabbatch);
                    end
                    movefile(out_file,paths.tmp_nii);
                end
            end
        end
    end
end





%% Motion correction
run_count = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        run_count = run_count + 1;
        %% Filename
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);
        if reslice_flags.unwarp && has_fmap(run_count)
            bold_file = fullfile(paths.tmp_nii,['unwarped_' bold_name ext]);
        end

        RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
        load(RealignRef_file)
        
        %% Perform realignment for this run
        out_base = fullfile(paths.tmp_nii,['MC_' bold_name]);
        out_file = [out_base '.nii'];
        par_out_file = [out_base '.par'];
        rp_file = fullfile(paths.rp_files,['rp_' bold_name '.txt']);
        
        if ~(exist(out_file,'file') == 2) || ~same_analysis
            fprintf('\nCorrecting motion for session %i, run %i\n',session,run)
            if ~exist(paths.rp_files,'dir');mkdir(paths.rp_files);end % create folder if non-existant
            if exist(out_file,'file') == 2; delete(out_file);end
            if strcmp(realign_FSL.MC_method,'mcflirt')
                system(sprintf('%smcflirt -in %s -out %s -refvol %i -plots -cost %s -dof %i',paths.FSL_prefix,bold_file,out_base,ref_volume-1,realign_FSL.mcflirt_cost,realign_FSL.mcflirt_dof));
                movefile(par_out_file,rp_file)
            elseif strcmp(realign_FSL.MC_method,'ants')
                pre_rp_file = fullfile(paths.tmp_nii,['rp_' bold_name '.par']); % rp file of the first realignment estimate
                copyfile(pre_rp_file,rp_file)
                % average selected volumes
                P = spm_vol(bold_file); % read scan
                Y = spm_read_vols(P);
                avg_Y = mean(Y(:,:,:,selected_vols),4);
                avg_out_file = fullfile(paths.tmp_nii,['avg_for_ants_' bold_name '.nii']);
                P = P(1);
                P.fname = avg_out_file;
                P.n = [1 1];
                spm_write_vol(P,avg_Y); % write it
                % antsMotionCorr
                str = [fullfile('$ANTSPATH','antsMotionCorr') ' -d 3 ',...
                      '-o [' out_base ',' out_base '.nii.gz] ',... 
                      '-m MI[' avg_out_file ',' bold_file ',1,1,Random,0.05] ',...
                      '-t Rigid[0.005] -i 20 -u 1 -e 1 -s 0 -f 1 ',...
                      '-m MI[' avg_out_file ',' bold_file ',1,1,Random,0.05] ',...
                      '-t Affine[0.005] -i 20 -u 1 -e 1 -s 0 -f 1 ',...
                      '-m MI[' avg_out_file ',' bold_file ',1,2] ',...
                      '-t SyN[0.15,3,0.5] -i 20 -u 1 -e 1 -s 0 -f 1 -n 10'];
                      
                system(str);
            else
                error('realign_FSL.MC_method ''%s'' does not exist!!',realign_FSL.MC_method)
            end
            system(sprintf('gunzip %s',out_file));  
        end


        %% Put to double datatype
        % [~,datatype_out] = system([paths.FSL_prefix 'fslval ' out_file ' datatype']);
        % datatype_out = str2double(datatype_out);
        
        % if datatype_out ~= 64
        %     fprintf('Converting session %i run % to ''double''\n',session,run)
        %     system(sprintf('gzip %s',out_file));
        %     system(sprintf('%sfslmaths %s %s -odt 64',paths.FSL_prefix,out_file,out_file)); % apply values to resampled template
        %     system(sprintf('gunzip %s',out_file));
        % end
        
        %% ref_scans
        ref_scans{run_count,1} = out_file; % ref_scans become the MC ones
    end
end




%% Extract each ref vol & compute mean/std
fprintf('Computing run averages...\n')
M_avg = nan(size(ref_scans,1),1);
MS_avg = nan(size(ref_scans,1),1);
avg_out_files = cell(size(ref_scans,1),1);
for i = 1:size(ref_scans,1)
    [~,epi_name] = fileparts(ref_scans{i,1});
    fprintf('%s\n',epi_name)
    
    % make average of this run
    [~,name] = fileparts(ref_scans_orig{i,1});
    RealignRef_file = fullfile(paths.realign,[name '_RealignRef.mat']); % get ref file
    load(RealignRef_file) % get ref file
    P = spm_vol(ref_scans{i,1}); % read scan
    Y = spm_read_vols(P);
    avg_Y = mean(Y(:,:,:,selected_vols),4);
    
    avg_out_files{i} = fullfile(paths.tmp_nii,sprintf('%s_avg.nii',epi_name));
    P = P(1);
    P.fname = avg_out_files{i};
    P.n = [1 1];
    spm_write_vol(P,avg_Y); % write it
    
    % Mean / std
    [~,M] = system(sprintf('%sfslstats %s -M',paths.FSL_prefix,avg_out_files{i}));
    [~,S] = system(sprintf('%sfslstats %s -S',paths.FSL_prefix,avg_out_files{i}));
    M_avg(i) = str2double(M);
    MS_avg(i) = str2double(M) / str2double(S);
end


[~,rv] = min(abs(MS_avg - nanmean(MS_avg))); % Reference run



Master.ref_vol = ref_scans{rv,1};
[~,name] = fileparts(Master.ref_vol);
Master.session = str2double(name(strfind(name,'_ses-')+5:strfind(name,'_ses-')+6));
Master.session_id = find([SR.session] == Master.session);
Master.run = str2double(name(strfind(name,'_run-')+5:strfind(name,'_run-')+6));
Master.run_id = find(SR(Master.session_id).runs == Master.run);
Master.vol = SR(Master.session_id).ref_volume(Master.run_id);

save(paths.reference_scan_infos,'Master')

fprintf(fid,'\n\nReference for all alignments:\nSession %i, run %i (volume %i)\n\n',Master.session,Master.run,Master.vol);
fclose(fid);

fprintf('Done\n\n')

%% Write reference volume of the whole analysis
if exist(paths.reference_scan,'file') == 2; delete(paths.reference_scan);end
copyfile(avg_out_files{rv},paths.reference_scan)





%% Setting paths & -p option for bash scripts
debias_path = which('T1xT2BiasFieldCorrection.sh');
BET_path = which('T1xT2BET.sh');
IterREGBET_path = which('IterREGBET.sh');
if isfield(paths,'FSL_prefix') && ~isempty(paths.FSL_prefix)
    popt = ['-p ' paths.FSL_prefix];
else
    popt ='';
end






%% Apply bias field correction to the reference volumes
debiased_files = cell(length(avg_out_files),1);
for i = 1:length(avg_out_files)
    [~,epi_name] = fileparts(avg_out_files{i});
    in_file = avg_out_files{i};
    debiased_files{i} = fullfile(paths.tmp_nii,sprintf('%s_debiased.nii.gz',epi_name));
    system(sprintf('bash %s -t1 %s -t2 %s -bet 2 -f 0.1 %s',debias_path,in_file,in_file,popt));
end



if realign_FSL.IR_BET
    %% BETing debiased averages
    BETed_files = cell(length(debiased_files),1);
    for i = 1:length(debiased_files)
        [~,epi_name,ext] = fileparts(debiased_files{i});
        if strcmp(ext,'.gz'); [~,epi_name] = fileparts(epi_name); end
        in_file = debiased_files{i};
        BETed_files{i} = fullfile(paths.tmp_nii,sprintf('%s_BET.nii.gz',epi_name));
        system(sprintf('bash %s -t1 %s -t2 %s -n 2 -f %f %s',BET_path,in_file,in_file,realign_FSL.IR_BETf,popt));
    end




    %% Apply anat brain mask to debiased averages (IterREGBET)
    brain_files = cell(length(debiased_files),1);
    for i = 1:length(debiased_files)
        [~,epi_name,ext] = fileparts(debiased_files{i});
        if strcmp(ext,'.gz'); [~,epi_name] = fileparts(epi_name); end
        inw = debiased_files{i};
        inb = BETed_files{i};
        refb = paths.anat.brain;
        brain_files{i} = fullfile(paths.tmp_nii,sprintf('%s_IRbrain.nii.gz',epi_name));
        system(sprintf('bash %s -inw %s -inb %s -refb %s -dof 6 %s',IterREGBET_path,inw,inb,refb,popt));
    end

    files_to_reg = brain_files;
else
    files_to_reg = debiased_files;
end






%% First registration of reference volumes
fprintf('Preliminary Inter run registration\n')
first_reg_brain_files = cell(length(files_to_reg),1);
first_reg_whole_files = cell(length(files_to_reg),1);
ref_brain_file = files_to_reg{rv};
ref_whole_file = debiased_files{rv};
for i = 1:length(files_to_reg)
    fprintf('Registering run %i / %i...\n',i,length(files_to_reg))
    [~,this_epi_name] = fileparts(ref_scans{i,1});
    in_brain_file = files_to_reg{i};
    in_whole_file = debiased_files{i};
    first_reg_brain_files{i} = fullfile(paths.tmp_nii,['first_reg_brain_' this_epi_name '.nii.gz']);
    first_reg_whole_files{i} = fullfile(paths.tmp_nii,['first_reg_whole_' this_epi_name '.nii.gz']);
    if i ~= rv
        if strcmp(realign_FSL.IR_method,'flirt')
            system(sprintf('%sflirt -in %s -ref %s -dof %i -cost %s -out %s',paths.FSL_prefix,in_brain_file,ref_brain_file,realign_FSL.flirt_dof,realign_FSL.flirt_cost,first_reg_brain_files{i})); % compute xfm
        elseif strcmp(realign_FSL.IR_method,'ants')
            ANTs_out_base = fullfile(paths.tmp_nii,['SyN_' this_epi_name '_']);
            if realign_FSL.IR_BET
                system(sprintf('%s -d 3 -f %s -f %s -m %s -m %s -o %s -t sr -j 1',fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh'),ref_brain_file,ref_whole_file,in_brain_file,in_whole_file,ANTs_out_base));
            else
                system(sprintf('%s -d 3 -f %s -m %s -o %s -t sr -j 1',fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh'),ref_brain_file,in_brain_file,ANTs_out_base));
            end
            copyfile([ANTs_out_base 'Warped.nii.gz'],first_reg_brain_files{i})
            warp_file = [ANTs_out_base '1Warp.nii.gz'];
            xfm_file = [ANTs_out_base '0GenericAffine.mat'];
            system(sprintf('%s -i %s -r %s -o %s -t %s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),in_whole_file,ref_whole_file,first_reg_whole_files{i},warp_file,xfm_file));
        else
            error('realign_FSL.IR_method ''%s'' does not exist!!',realign_FSL.IR_method)
        end
    else
        copyfile(in_brain_file,first_reg_brain_files{i})
        copyfile(in_whole_file,first_reg_whole_files{i})
    end

    fprintf('Done\n')
end





%% Merge ref volumes and average them
fprintf('Merging and averaging ref brain volumes...')
merge_in_files = sprintf('%s ',first_reg_brain_files{:});
merge_out_file = fullfile(paths.resliced,'merged_ref_brain_vols.nii.gz');
system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,merge_out_file,merge_in_files)); % merge anat runs
ref_brain_file = fullfile(paths.resliced,'Inter_run_registration_brain_reference.nii.gz');
system(sprintf('%sfslmaths %s -Tmean %s',paths.FSL_prefix,merge_out_file,ref_brain_file)); % average merged anat runs
fprintf('Done\n\n')

if strcmp(realign_FSL.IR_method,'ants')
    fprintf('Merging and averaging ref whole volumes...')
    merge_in_files = sprintf('%s ',first_reg_whole_files{:});
    merge_out_file = fullfile(paths.resliced,'merged_ref_whole_vols.nii.gz');
    system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,merge_out_file,merge_in_files)); % merge anat runs
    ref_whole_file = fullfile(paths.resliced,'Inter_run_registration_whole_reference.nii.gz');
    system(sprintf('%sfslmaths %s -Tmean %s',paths.FSL_prefix,merge_out_file,ref_whole_file)); % average merged anat runs
    fprintf('Done\n\n')
end




%% Final registration of reference volumes
fprintf('Final Inter run registration\n')
for i = 1:length(files_to_reg)
    fprintf('Registering run %i / %i...\n',i,length(files_to_reg))
    [~,this_epi_name] = fileparts(ref_scans{i,1});
    in_brain_file = files_to_reg{i};
    in_whole_file = debiased_files{i};
    in_file_4D = ref_scans{i,1};
    [~,epi_name,ext] = fileparts(EPIs{i,1});
    out_file_4D = fullfile(paths.resliced,[reslice_flags.prefix epi_name ext]);
    if exist(out_file_4D,'file') == 2; delete(out_file_4D);end
    if strcmp(realign_FSL.IR_method,'flirt')
        omat = fullfile(paths.tmp_nii,sprintf('%i_to_%i.xfm',i,rv));
        system(sprintf('%sflirt -in %s -ref %s -omat %s -dof %i -cost %s',paths.FSL_prefix,in_brain_file,ref_brain_file,omat,realign_FSL.flirt_dof,realign_FSL.flirt_cost)); % compute xfm
        fprintf('Applying transform to time series...')
        system(sprintf('%sflirt -in %s -ref %s -init %s -applyxfm -out %s',paths.FSL_prefix,in_file_4D,ref_brain_file,omat,out_file_4D)); % apply to time series
    elseif strcmp(realign_FSL.IR_method,'ants')
        ANTs_out_base = fullfile(paths.tmp_nii,['SyN_' this_epi_name '_']);
        if realign_FSL.IR_BET
            system(sprintf('%s -d 3 -f %s -f %s -m %s -m %s -o %s -t sr -j 1',fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh'),ref_brain_file,ref_whole_file,in_brain_file,in_whole_file,ANTs_out_base));
        else
            system(sprintf('%s -d 3 -f %s -m %s -o %s -t sr -j 1',fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh'),ref_brain_file,in_brain_file,ANTs_out_base));
        end
        warp_file = [ANTs_out_base '1Warp.nii.gz'];
        xfm_file = [ANTs_out_base '0GenericAffine.mat'];
        fprintf('Applying warp to time series...')
        system(sprintf('%s -i %s -e 3 -r %s -o %s -t %s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),in_file_4D,ref_whole_file,out_file_4D,warp_file,xfm_file));
    end
    system(sprintf('gunzip %s',out_file_4D));
    fprintf('Done\n')
end



%% Construct a verification scan & average
% The verification scan comprises 3 steady volumes of each run
fprintf('\nConstructing verification scan\n')

if exist('avg_Y','var');clear avg_Y;end

verif_vol = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        fprintf('Session %i, run %i',session,run)
        bold_file = SR(s).filename{r};
        [~,name,ext] = fileparts(bold_file);
        % get ref vol
        RealignRef_file = fullfile(paths.realign,[name '_RealignRef.mat']);
        load(RealignRef_file)
        
        if any(selected_vols)
            % read scan
            resliced_scan = fullfile(paths.resliced,[reslice_flags.prefix name ext]);
            P = spm_vol(resliced_scan);
            Y = spm_read_vols(P);
            % select verification images
            a = 1:numel(P);
            ok_vols = a(selected_vols);
            v = ok_vols(round(linspace(1,length(ok_vols),mvt_params.n_vols_verif))); % select n_vols_verif volumes to put in verif scan
            % Average
            if ~exist('grand_avg_Y','var')
                grand_avg_Y = mean(Y(:,:,:,ok_vols),4);
            else
                grand_avg_Y = grand_avg_Y + mean(Y(:,:,:,ok_vols),4);
            end
            % write vols
            nvols = length(v);
            filename = fullfile(paths.resliced,sprintf('Verif_%s.nii',paths.realign_name));
            if (s == 1) && (r == 1)
                verifmat = P(1).mat;
            end
            for i = 1:nvols
                fprintf('.')
                verif_vol = verif_vol +1;
                VO = P(v(i));
                VO.fname = filename;
                VO.n = [verif_vol 1];
                VO.mat = verifmat;
                % add session and run number to each volume
                txtstr = sprintf('%i - %i',session,run);
                txtval = max(Y(:));
                vol2write = zeros(size(Y,1),size(Y,2),size(Y,3));
                for z = 1:size(Y,3)
                    slice = Y(:,:,z,v(i));
                    vol2write(:,:,z) = addtext2slice(slice,txtstr,txtval);
                end
                % write it
                spm_write_vol(VO,vol2write);
            end
        end
        fprintf('\n')
    end
end

%% Write average
fprintf('\nWriting average scan\n')
grand_avg_Y = grand_avg_Y / n_total_runs;
VO.fname = paths.average;
VO.n = [1 1];
spm_write_vol(VO,grand_avg_Y);


%% Inter session cost
% inter_session_cost


%% Remove tmp_nii folder
if exist(paths.tmp_nii,'dir')
    fprintf('\nCleaning temp folder\n')
    delete(fullfile(paths.tmp_nii,'*'))
    rmdir(paths.tmp_nii)
end














%% TRASH

% To place around line 26 if you want to T1 normalize before estimating
% realignment:

%         %% Check whether the TR is constant or variable. Create nomalized series if TR is variable
%         if (max(unique(diff(scan_log(s,r).Scan_trigs))) / scan_info(s,r).RepetitionTime) > 1.1
%             T1norm_file = fullfile(paths.T1_normalized,['T1norm_' name ext]);
%
%             if ~exist(T1norm_file,'file')
%                 fprintf('T1-normalizing session %i, run %i...\n',session,run)
%                 P = spm_vol(bold_file);
%                 Y = spm_read_vols(P);
%
%                 % load PCA volumes
%                 no_mask_file = fullfile(paths.tmp_nii,[name '_no-mask' ext]);
%                 system(sprintf('%sfslroi %s %s 0 1',paths.FSL_prefix,bold_file,no_mask_file));
%                 system(sprintf('%sfslmaths %s -bin %s',paths.FSL_prefix,no_mask_file,no_mask_file));
%                 system(sprintf('gunzip %s',no_mask_file));
%                 [PCA_mat_file,PCA_nii_file] = select_PCA_file(bold_file,paths,'no_mask');
%                 if ~exist(PCA_mat_file,'file')
%                     run_PCA_image(bold_file,paths,'no_mask',no_mask_file);
%                 end
%
%                 fmristat_PCA = load(PCA_mat_file); % load fMRISTAT PCA components
%
%                 Ppca = spm_vol(PCA_nii_file);
%                 Ypca = spm_read_vols(Ppca);
%
%                 PC = 1; % only consider the first PC
%                 Ypca1 = squeeze(Ypca(:,:,:,PC));
%
%                 % find the RMS ratio between min and max RMS for each volume of our time series Y
%                 maxRMS = 0;
%                 minRMS = inf;
%                 for t = 1:size(Y,4)
%                     yt = squeeze(Y(:,:,:,t));
%                     RMS = sqrt(mean(yt(Ypca1 > 0.99).^2));
%                     if RMS > maxRMS
%                         maxRMS = RMS;
%                     end
%                     if RMS < minRMS
%                         minRMS = RMS;
%                     end
%                 end
%                 RMSratio = minRMS / maxRMS;
%
%                 % Normalize component with RMSratio (the output w is the factor to apply to the time series in order to temporally normalize it)
%                 v = -fmristat_PCA.PCs(:,PC);
%                 v = v + (1 - max(v));
%                 minv = min(v);
%                 w = interp1([minv 1],[RMSratio 1],v);
%
%                 Ypca(isnan(Ypca)) = 0; % convert NaNs to zeros
%
%                 % Apply temporal + spatial normalization (using PCA image)
%                 for t = 1:size(Y,4)
%                     Ypcaw = Ypca(:,:,:,PC);
%                     Ypcaw = interp1([1 0 -1],[w(t) 1 (1/w(t))],Ypcaw); % weight PCA image with w
% %                     Y(:,:,:,t) = Y(:,:,:,t) .* Ypcaw; % normalize (both spatial & temporal)
%                     Y(:,:,:,t) = Y(:,:,:,t) .* w(t); % no spatial normalization (temporal only)
%                 end
%                 clear Ypcaw Ypca
%
%                 fprintf('Writing T1-normalized series...\n')
%                 new_P = P;
%                 if ~exist(paths.T1_normalized,'dir');mkdir(paths.T1_normalized);end % create folder if non-existant
%                 for i = 1:size(Y,4)
%                     new_P(i) = P(1);
%                     new_P(i).fname = T1norm_file;
%                     new_P(i).n = [i 1];
%                     spm_write_vol(new_P(i),squeeze(Y(:,:,:,i)));
%                 end
%             end
%             bold_file = T1norm_file;
%         end
