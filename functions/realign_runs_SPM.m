%% PARAMETERS

%%%%%% Run 'set_parameters' beforehand

%% SCRIPT
if ~exist(paths.tmp_nii,'dir');mkdir(paths.tmp_nii);end % create tmp folder if non-existant
fid = fopen(paths.log_file_realign,'w');

%% Check if this particular realignment already exists
same_analysis = 0;
if exist(paths.realign_params,'file') == 2
    load(paths.realign_params)
    if isequaln(realign_flags,saved_realign_flags) && isequaln(mvt_params,saved_mvt_params)
        same_analysis = 1;
    end
end

%% Find representative scans of each run
M = struct('mat',[]);
ref_scans = cell(n_total_runs,1);
run_count = 0;
vols_selected = zeros(n_total_runs,5);
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    SR(s).ref_volume = nan(length(runs),1);
    SR(s).selected_vols = cell(length(runs),1);
    SR(s).tmp_filename = cell(length(runs),1);
    
    for r = 1:length(runs)
        run = runs(r);
        %% Filename
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);
        if strcmp(ext,'.gz'); [~,bold_name,ext] = fileparts(bold_name); end
        
        %% Copy original bold files to the tmp folder
        fprintf('Copying files to tmp folder...')
        tmp_bold_file = fullfile(paths.tmp_nii,[bold_name ext]);
        copyfile(bold_file,tmp_bold_file);
        SR(s).tmp_filename{r} = tmp_bold_file;
        fprintf(' Done\n')
        
        %% List scans
        Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
        for vol = 1:numel(Scans)
            Scans{vol,1} = sprintf('%s,%i',bold_file,vol);
        end
        Scans = char(Scans); % without this line, spm thinks you have nvolumes sessions
        
        %% Check if this particular realignment has already been saved
        RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
        if ~(exist(RealignRef_file,'file') == 2) || ~same_analysis % if the file does not exist, perform realignment estimate
            %% Find ref vol
            fprintf('\nFinding representative volume for session %i, run %i:\n',session,run)
            [ref_volume, selected_vols, Q, i1, i2] = find_refvol_SPM(Scans,realign_flags,mvt_params,paths);
            
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
            save(RealignRef_file,'ref_volume','selected_vols')
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
        ref_scans{run_count} = Scans(ref_volume,:);
        
        %% Save vols selected
        if exist('vols_selected','var'); clear vols_selected; end
        vols_selected(run_count,1) = session;
        vols_selected(run_count,2) = run;
        vols_selected(run_count,3) = sum(selected_vols);
        vols_selected(run_count,4) = length(selected_vols) - sum(selected_vols);
        vols_selected(run_count,5) = length(selected_vols);
    end
end

%% Save vols selected
vols_sel_file = fullfile(paths.realign,['sub-' paths.subject '_task-' paths.task '_VolsSelected']);
save(vols_sel_file,'vols_selected')

saved_realign_flags = realign_flags;
saved_mvt_params = mvt_params;
save(paths.realign_params,'saved_realign_flags','saved_mvt_params')


%% Realign & Unwarp (Calculate VDM if fmap present)
spm('defaults', 'FMRI');
spm_jobman('initcfg'); % initialization
spm_get_defaults('cmdline',true)
if reslice_flags.unwarp
    has_fmap = zeros(run_count,1);
    all_VDMs = cell(run_count,1);
    run_count_fm = 0;
    for s = 1:numel(SR)
        session = SR(s).session;
        runs = SR(s).runs;
        if exist('fmap','var'); clear fmap; end
        fieldmaps = spm_BIDS(BIDS,'data','sub',paths.subject,'modality','fmap','type','fieldmap','ses',sprintf('%02.0f',session));
        if ~isempty(fieldmaps)
            SR(s).fieldmap = cell(length(runs),1);
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
                        if ismember(epi_run,runs) % if current run is part of this analysis

                            bold_files = spm_BIDS(BIDS,'data','sub',paths.subject,'ses',sprintf('%02.0f',session),'task',paths.task,'type','bold');
                            epi_pathstr = fileparts(bold_files{1});
                            epi_json_file = fullfile(epi_pathstr,[epi_name '.json']);
                            if exist(epi_json_file,'file')
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
                run_count_fm = run_count_fm + 1;
                found_run = 0;
                for ci = 1:numel(IntendedFor)
                    c = ~cellfun(@isempty,strfind(IntendedFor{ci},sprintf('sub-%s_ses-%02.0f_task-%s_run-%02.0f',paths.subject,session,paths.task,run)));
                    if any(c)
                        found_run = 1;
                        break
                    end
                end
                if found_run
                    has_fmap(run_count_fm) = 1;
                    all_VDMs{run_count_fm} = SR(s).vdmfile{r};
                end
            end
        end
    end
    
    
    if sum(has_fmap) >= (length(has_fmap) / 2) % there is a majority of runs with a fieldmap
        candidate_vols = ref_scans(logical(has_fmap));
        master_has_fmap = 1;
        VDMs = all_VDMs(logical(has_fmap));
    else % there is a majority of runs without fieldmap
        candidate_vols = ref_scans(logical(~has_fmap));
        master_has_fmap = 0;
    end
else % if reslice_flags.unwarp
    candidate_vols = ref_scans;
end



%% Choose master reference volume based on fieldmaps presence
P = spm_realign(char(candidate_vols),realign_flags);
n = length(P);
Q = zeros(n,6);
for j=1:n
    qq     = spm_imatrix(P(j).mat/P(1).mat);
    Q(j,:) = qq(1:6);
end
Q(:,4:6) = rad2deg(Q(:,4:6)); % radians to degrees
rv = find_mean_vol(Q);
Master.ref_vol = candidate_vols{rv};
[~,bold_name] = fileparts(Master.ref_vol);
Master.session = str2double(bold_name(strfind(bold_name,'_ses-')+5:strfind(bold_name,'_ses-')+6));
Master.session_id = find([SR.session] == Master.session);
Master.run = str2double(bold_name(strfind(bold_name,'_run-')+5:strfind(bold_name,'_run-')+6));
Master.run_id = find(SR(Master.session_id).runs == Master.run);
Master.vol = SR(Master.session_id).ref_volume(Master.run_id);

save(paths.reference_scan_infos,'Master')

fprintf(fid,'\n\nReference for all alignments:\nSession %i, run %i, volume %i\n\n',Master.session,Master.run,Master.vol);
fclose(fid);

%% Generate figure
figureprep([0 0 1280 960]);
realign_plot_Q(Q,'Master volume identification')
hold on
plot([rv rv],[min(min(Q(:,4:6))) max(max(Q(:,4:6)))],'LineWidth',2)
text(rv+1,(max(max(Q(:,4:6))) / 2),sprintf('Master\nvolume:\nsession %i\nrun %i',Master.session,Master.run))
hold off
fig_prefix = fullfile(paths.realign,'Master_volume_identificatiion');
figurewrite(fig_prefix,[],0,paths.realign); % the 0 is to force eps figure

%% Write reference volume of the whole analysis
if ~exist(paths.resliced,'dir')
    mkdir(paths.resliced) % create folder if non-existant
else
    delete(fullfile(paths.resliced,'*'))
end

system(sprintf('%sfslroi %s %s %i 1',paths.FSL_prefix,SR(Master.session_id).filename{Master.run_id},paths.reference_scan,Master.vol-1)); % fslroi indexing starts at 0, hence the '-1' here
system(sprintf('gunzip %s',paths.reference_scan));




if reslice_flags.unwarp && any(has_fmap)
    %% Realign & Unwarp
    % We have to do two separate realign & unwarp for sessions that have
    % fieldmaps and sessions that don't. From SPM manual:
    % Only add similar session data to a realign+unwarp branch, i.e., choose Data or Data+phase map
    % for all sessions, but don?t use them interchangeably.
    if exist('matlabbatch','var'); clear matlabbatch; end
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = realign_flags.quality;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = realign_flags.sep;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = realign_flags.fwhm;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = realign_flags.interp;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = realign_flags.wrap;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = realign_flags.fwhm;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [reslice_flags.which reslice_flags.mean];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = reslice_flags.interp;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = reslice_flags.wrap;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = reslice_flags.mask;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = reslice_flags.prefix;
    
    % Sessions that have fieldmaps
    if master_has_fmap % Master reference scan has an associated fieldmap > at least the majority of the runs also have an associated fielmap
        master_found = 0;
        for d = 1:sum(has_fmap) % first treat runs that have an associated fieldmap
            if d == rv
                i = 1; % place master reference run at the beginning
                master_found = 1;
            else
                if ~master_found
                    i = d + 1;
                else
                    i = d;
                end
            end
            scan_file = candidate_vols{d};
            pmscan_file = VDMs{d};
            session = str2double(scan_file(strfind(scan_file,'_ses-')+5:strfind(scan_file,'_ses-')+6));
            s = find([SR.session] == session);
            run = str2double(scan_file(strfind(scan_file,'_run-')+5:strfind(scan_file,'_run-')+6));
            r = find(SR(s).runs == run);
            
            Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
            for vol = 1:numel(Scans)
                Scans{vol,1} = sprintf('%s,%i',SR(s).tmp_filename{r},vol);
            end
            
            Scans = [Scans(SR(s).ref_volume(r),:);Scans]; % put the master ref_vol (ref_vol of 1st run / 1st session) at the beginning of the list. Every reslice operations will be done relative to this master ref_vol
            Scans(SR(s).ref_volume(r) + 1) = []; % keep ref vol at the beginning of the list but remove it in the middle
            
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans = Scans;
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = {pmscan_file};
        end
        
        % Run job
        output_list = spm_jobman('run', matlabbatch);
        
        if any(~has_fmap) % if some runs don't have associated fieldmap
            no_fmap_vols = ref_scans(logical(~has_fmap));
            % Copy reference scan to tmp folder
            copyfile(paths.reference_scan,fullfile(paths.tmp_nii,'Reference_scan.nii'));
            % Place the copy of reference scan at the first session
            matlabbatch{1}.spm.spatial.realignunwarp.data = [];
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans = {fullfile(paths.tmp_nii,'Reference_scan.nii,1')};
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = '';
            for d = 1:sum(~has_fmap)
                scan_file = no_fmap_vols{d};
                session = str2double(scan_file(strfind(scan_file,'_ses-')+5:strfind(scan_file,'_ses-')+6));
                s = find([SR.session] == session);
                run = str2double(scan_file(strfind(scan_file,'_run-')+5:strfind(scan_file,'_run-')+6));
                r = find(SR(s).runs == run);
                
                Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
                for vol = 1:numel(Scans)
                    Scans{vol,1} = sprintf('%s,%i',SR(s).tmp_filename{r},vol);
                end
                
                Scans = [Scans(SR(s).ref_volume(r),:);Scans]; % put the master ref_vol (ref_vol of 1st run / 1st session) at the beginning of the list. Every reslice operations will be done relative to this master ref_vol
                Scans(SR(s).ref_volume(r) + 1) = []; % keep ref vol at the beginning of the list but remove it in the middle
                
                matlabbatch{1}.spm.spatial.realignunwarp.data(d+1).scans = Scans;
                matlabbatch{1}.spm.spatial.realignunwarp.data(d+1).pmscan = '';
            end
            % Run job
            output_list = spm_jobman('run', matlabbatch);
        end
    else % Master reference scan does not have an associated fieldmap > at least the majority of the runs also don't have an associated fielmap
        master_found = 0;
        for d = 1:sum(~has_fmap) % first treat runs that don't have an associated fieldmap
            if d == rv
                i = 1; % place master reference run at the beginning
                master_found = 1;
            else
                if ~master_found
                    i = d + 1;
                else
                    i = d;
                end
            end
            scan_file = candidate_vols{d};
            session = str2double(scan_file(strfind(scan_file,'_ses-')+5:strfind(scan_file,'_ses-')+6));
            s = find([SR.session] == session);
            run = str2double(scan_file(strfind(scan_file,'_run-')+5:strfind(scan_file,'_run-')+6));
            r = find(SR(s).runs == run);
            
            Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
            for vol = 1:numel(Scans)
                Scans{vol,1} = sprintf('%s,%i',SR(s).tmp_filename{r},vol);
            end
            
            Scans = [Scans(SR(s).ref_volume(r),:);Scans]; % put the master ref_vol (ref_vol of 1st run / 1st session) at the beginning of the list. Every reslice operations will be done relative to this master ref_vol
            Scans(SR(s).ref_volume(r) + 1) = []; % keep ref vol at the beginning of the list but remove it in the middle
            
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans = Scans;
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = '';
        end
        
        % Run job
        output_list = spm_jobman('run', matlabbatch);
        
        if any(has_fmap) % if some runs have an associated fieldmap
            fmap_vols = ref_scans(logical(has_fmap));
            % Copy reference scan to tmp folder
            copyfile(paths.reference_scan,fullfile(paths.tmp_nii,'Reference_scan.nii'));
            % create zero vdm
            out_file = fullfile(paths.tmp_nii,'VDM_zeros.nii');
            if exist(out_file,'file') == 2; delete(out_file);end
            system(sprintf('%sfslmaths %s -mul 0 %s',paths.FSL_prefix,paths.reference_scan,out_file));
            system(sprintf('gunzip %s',out_file));
            % Place the copy of reference scan at the first session
            matlabbatch{1}.spm.spatial.realignunwarp.data = [];
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans = {fullfile(paths.tmp_nii,'Reference_scan.nii,1')};
            matlabbatch{1}.spm.spatial.realignunwarp.data(1).pmscan = {out_file};
            for d = 1:sum(has_fmap)
                scan_file = fmap_vols{d};
                pmscan_file = VDMs{d};
                session = str2double(scan_file(strfind(scan_file,'_ses-')+5:strfind(scan_file,'_ses-')+6));
                s = find([SR.session] == session);
                run = str2double(scan_file(strfind(scan_file,'_run-')+5:strfind(scan_file,'_run-')+6));
                r = find(SR(s).runs == run);
                
                Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
                for vol = 1:numel(Scans)
                    Scans{vol,1} = sprintf('%s,%i',SR(s).tmp_filename{r},vol);
                end
                
                Scans = [Scans(SR(s).ref_volume(r),:);Scans]; % put the master ref_vol (ref_vol of 1st run / 1st session) at the beginning of the list. Every reslice operations will be done relative to this master ref_vol
                Scans(SR(s).ref_volume(r) + 1) = []; % keep ref vol at the beginning of the list but remove it in the middle
                
                matlabbatch{1}.spm.spatial.realignunwarp.data(d+1).scans = Scans;
                matlabbatch{1}.spm.spatial.realignunwarp.data(d+1).pmscan = {pmscan_file};
            end
            % Run job
            output_list = spm_jobman('run', matlabbatch);
        end
    end
    
else % if reslice_flags.unwarp && any(has_fmap)
    %% Estimate & Reslice
    if exist('matlabbatch','var'); clear matlabbatch; end
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = realign_flags.quality;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = realign_flags.sep;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = realign_flags.fwhm;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = realign_flags.interp;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = realign_flags.wrap;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [reslice_flags.which reslice_flags.mean];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = reslice_flags.interp;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = reslice_flags.wrap;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = reslice_flags.mask;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = reslice_flags.prefix;
    
    master_found = 0;
    for d = 1:run_count
        if d == rv
            i = 1; % place master reference run at the beginning
            master_found = 1;
        else
            if ~master_found
                i = d + 1;
            else
                i = d;
            end
        end
        scan_file = candidate_vols{d};
        session = str2double(scan_file(strfind(scan_file,'_ses-')+5:strfind(scan_file,'_ses-')+6));
        s = find([SR.session] == session);
        run = str2double(scan_file(strfind(scan_file,'_run-')+5:strfind(scan_file,'_run-')+6));
        r = find(SR(s).runs == run);
        
        Scans = cell(scan_info(s,r).NumberOfVolumesInFile,1);
        for vol = 1:numel(Scans)
            Scans{vol,1} = sprintf('%s,%i',SR(s).tmp_filename{r},vol);
        end
        
        Scans = [Scans(SR(s).ref_volume(r),:);Scans]; % put the master ref_vol (ref_vol of 1st run / 1st session) at the beginning of the list. Every reslice operations will be done relative to this master ref_vol
        Scans(SR(s).ref_volume(r) + 1) = []; % keep ref vol at the beginning of the list but remove it in the middle
        
        matlabbatch{1}.spm.spatial.realign.estwrite.data{i} = Scans;
    end
    
    % Run job
    output_list = spm_jobman('run', matlabbatch);
end



%% Move resliced files
source_to_move = fullfile(paths.tmp_nii,sprintf('%ssub-%s*.nii',reslice_flags.prefix,paths.subject));
movefile(source_to_move,paths.resliced); % move resliced scans to the resliced directory

source_to_move = fullfile(paths.tmp_nii,sprintf('rp_sub-%s*.txt',paths.subject));
if ~exist(paths.rp_files,'dir');mkdir(paths.rp_files);end % create folder if non-existant
movefile(source_to_move,paths.rp_files); % move resliced scans to the rp_files directory



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
        [~,bold_name,ext] = fileparts(bold_file);
        % get ref vol
        RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
        load(RealignRef_file)
        
        if any(selected_vols)
            % read scan
            resliced_scan = fullfile(paths.resliced,[reslice_flags.prefix bold_name ext]);
            P = spm_vol(resliced_scan);
            Y = spm_read_vols(P);
            % select verification images
            a = 1:numel(P);
            ok_vols = a(selected_vols);
            v = ok_vols(round(linspace(1,length(ok_vols),mvt_params.n_vols_verif))); % select n_vols_verif volumes to put in verif scan
            % Average
            if ~exist('avg_Y','var')
                avg_Y = mean(Y(:,:,:,ok_vols),4);
            else
                avg_Y = avg_Y + mean(Y(:,:,:,ok_vols),4);
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

clear Y P

%% Write average
avg_Y = avg_Y / n_total_runs;
VO.fname = paths.average;
VO.n = [1 1];
spm_write_vol(VO,avg_Y);

%% Remove tmp_nii folder
if exist(paths.tmp_nii,'dir')
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
