function anat_files = get_anat(BIDS,paths,anat_suffix,anat_session)

anat_files = fullfile(paths.anat_dir,sprintf('sub-%s_ses-%02.0f_%s_anat-files.json',paths.subject,anat_session,anat_suffix));

if ~exist(anat_files,'file')
    anat_file = spm_BIDS(BIDS,'data','sub',paths.subject,'type',anat_suffix,'ses',sprintf('%02.0f',anat_session));
    if ~exist(paths.anat_dir,'dir');mkdir(paths.anat_dir);end % create folder if non-existant

    if numel(anat_file) == 1
        anat_json.raw = anat_file{1};
    elseif isempty(anat_file)
        error('No anatomic scan of type''%s'' found in session %02.0f.\nModify your parameters file.',anat_suffix,anat_session)
    else % several anat files
        anat_runs = spm_BIDS(BIDS,'runs','sub',paths.subject,'type',anat_suffix,'ses',sprintf('%02.0f',anat_session),'modality','anat');
        if numel(anat_runs) > 1 % if several anat runs, average them
            fprintf('\nAveraging %s runs...\n',anat_suffix)
            anat_file = cell(numel(anat_runs),1);
            for r = 1:numel(anat_runs)
                a = spm_BIDS(BIDS,'data','sub',paths.subject,'type',anat_suffix,'ses',sprintf('%02.0f',anat_session),'modality','anat','run',anat_runs{r});
                anat_file{r} = a{1};
            end
            [~,a_name,ext] = fileparts(anat_file{1});
            if strcmp(ext,'.gz'); [~,a_name] = fileparts(a_name); end
            strtofind = ['_run-' anat_runs{1}];
            idx = strfind(a_name,strtofind);
            a_name(idx:idx+length(strtofind)-1) = '';
            mean_out_file = fullfile(paths.anat_dir,[a_name '_Mean.nii.gz']);
            
            merge_in_files = sprintf('%s ',anat_file{:});
            system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,mean_out_file,merge_in_files)); % merge anat runs
            system(sprintf('%sfslmaths %s -Tmean %s',paths.FSL_prefix,mean_out_file,mean_out_file)); % average merged anat runs
            
            % [json,step] = get_preproc_history(anat_file{1});
            % pstep = sprintf('PV_preprocessing_step%02.0f',step + 1);
            % json.(pstep).type = 'averaging';
            % json.(pstep).method = 'fsl';
            % json.(pstep).input = anat_file{1};
            % json.(pstep).params.inputs = merge_in_files;
            
            % spm_jsonwrite([mean_out_file '.json'],json,struct('indent','    '));
            
            anat_json.raw = mean_out_file;
        else
            anat_json.raw = anat_file{1};
        end
    end
    spm_jsonwrite(anat_files,anat_json,struct('indent','    '));
end

anat_json = spm_jsonread(anat_files);

if ~isfield(anat_json,'name')
    [~,a_name,ext] = fileparts(anat_json.raw);
    if strcmp(ext,'.gz'); [~,a_name] = fileparts(a_name); end
    anat_json.name = a_name;
    spm_jsonwrite(anat_files,anat_json,struct('indent','    '));
end

