%% SCRIPT
% to be run after parameters_*

%% Run spm_BIDS
fprintf('Retrieving BIDS informations...')
BIDS = spm_BIDS(paths.dataset);
fprintf(' done.\n\n')

%% Check participants.tsv
if isempty(BIDS.participants)
    error('There must be a ''participants.tsv'' file in the BIDS folder!!')
elseif ~isfield(BIDS.participants,'realign_name')
    error('There must be a ''realign_name'' column in the ''participants.tsv'' file!!')
elseif ~isfield(BIDS.participants,'results_name')
    error('There must be a ''results_name'' column in the ''participants.tsv'' file!!')
elseif ~isfield(BIDS.participants,'results_sub_folder')
    error('There must be a ''results_sub_folder'' column in the ''participants.tsv'' file!!')    
else
    %% Get necessary infos & paths
    templates
    paths.results = cell(numel(paths.subjects),1);
    paths.results_name = cell(numel(paths.subjects),1);
    paths.results_multi = fullfile(paths.dataset,['second_level_analysis_' paths.analysis_name]);
    paths.contrast_agent = cell(numel(paths.subjects),1);
    fprintf('ID\tName\trealign_name\tresults_name\tresults_sub_folder\n')
    for p = 1:numel(paths.subjects)
        pn = find(strcmp(BIDS.participants.participant_id,['sub-' paths.subjects{p}]));
        if isempty(pn)
            error('Could not find participant in the ''participants.tsv'' file!!')
        elseif length(pn) > 1
            error('There are more than one participant with the ID ''%s'' in the ''participants.tsv'' file!!',['sub-' paths.subject])
        else
            species = BIDS.participants.species{pn};
            if p == 1
                previous_species = species;
            else
                if ~strcmp(previous_species,species)
                    error('You are attempting to combine subjects of different species!!')
                end
            end

            paths.contrast_agent{p} = BIDS.participants.contrast_agent{pn};
            
            fprintf('%s\t%s\t%s\t%s\t%s\n',...
                BIDS.participants.participant_id{pn},...
                BIDS.participants.name{pn},...
                BIDS.participants.realign_name{pn},...
                BIDS.participants.results_name{pn},...
                BIDS.participants.results_sub_folder{pn})

            paths_realign_name = BIDS.participants.realign_name{pn};
            paths.results_name{p} = BIDS.participants.results_name{pn};
            paths_results_sub_folder = BIDS.participants.results_sub_folder{pn};

            paths_subject_analysis = fullfile(paths.dataset,['analysis_sub-' paths.subjects{p}]);
            paths_analysis = fullfile(paths_subject_analysis,paths_realign_name);
            paths_results = fullfile(paths_analysis,['results_' paths.results_name{p}]);
            if strcmpi(paths_results_sub_folder,'none')
                paths.results{p} = paths_results;
            else
                paths.results{p} = fullfile(paths_results,paths_results_sub_folder);
            end
        end
    end
    paths.contrast_agent = unique(paths.contrast_agent);
    if numel(paths.contrast_agent) > 1
        error('Trying to combine subjects with different contrats agents!!')
    end

    %% Get contrast_names if absent
    if ~exist('contrast_names','var')
        stim_level = 1;
        fprintf('\nAnalysis will be performed at the condition level for L%i\n',contrasts.L)
        % create sessions / runs structure
        SR = struct_sess_run(BIDS,'all',1,paths.subjects{1},paths.task);
        s = numel(SR); % get infos from last session
        session = SR(s).session;
        runs = SR(s).runs;
        r = 1; % get infos from first run
        run = runs(r);
        ScanLog_file = fullfile(paths.dataset,'sourcedata',['sub-' paths.subjects{1}],sprintf('ses-%02.0f',session),'func',[SR(s).namebase{r} '_ScanLog.mat']);
        if exist('ScanLog','var'); clear ScanLog; end
        load(ScanLog_file)
        L = contrasts.L;
        % L = numel(fieldnames(ScanLog.fMRISTAT)); % stim level is always last level
        contrasts.nconds = numel(ScanLog.fMRISTAT.(['L' num2str(L)]).names) - 1; % exclude silence
        contrast_names = cell(contrasts.nconds,1);
        for c = 1:contrasts.nconds % each condition vs silence
            contrast_names{c} = [ScanLog.fMRISTAT.(['L' num2str(L)]).names{c} '_vs_' ScanLog.fMRISTAT.(['L' num2str(L)]).names{end}];
        end
        ScanLog_ref = ScanLog;
        fprintf('\nYou can now launch ''run_second_level'' or ''brain_RDMs''\n')
    else
        stim_level = 0;
        fprintf('\nAnalysis will be performed for the following contrasts:\n')
        fprintf('%s\n',contrast_names{:})
        fprintf('\nYou can now launch ''run_second_level''\n')
    end
end










