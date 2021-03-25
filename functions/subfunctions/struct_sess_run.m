function SR = struct_sess_run(BIDS,sessions,runs,subject,task,scan_type,scan_modality)
% Construct the sessions/runs SR structure

if ~exist('scan_type','var')
    scan_type = 'bold';
end

if ~exist('scan_modality','var')
    scan_modality = 'func';
end

% Check if there are data
if isempty(task)
    C = spm_BIDS(BIDS,'data','sub',subject,'type',scan_type,'modality',scan_modality);
else
    C = spm_BIDS(BIDS,'data','sub',subject,'type',scan_type,'modality',scan_modality,'task',task);
end
if isempty(C)
    SR = [];
    return
end

% Get sessions from BIDS if sessions = 'all'
if ischar(sessions)
    if strcmp(sessions,'all')
%         sessions = spm_BIDS(BIDS,'sessions','sub',subject,'type',scan_type,'task',task); % query for 'sessions' does not work....... Using homemade function below instead
        if isempty(task)
            sessions = BIDS_sessions_query(BIDS,subject,scan_modality,scan_type);
        else
            sessions = BIDS_sessions_query(BIDS,subject,scan_modality,scan_type,task);
        end
        sessions = str2double(sessions);
    else
        error('In parameters: if ''sessions'' is a char array it must be sessions = ''all''')
    end
end

% Create SR variable
if iscell(runs)
    for s = 1:length(sessions)
        SR(s).session = sessions(s);
        SR(s).runs = runs{s}(:);
    end
elseif isnumeric(runs)
    for s = 1:length(sessions)
        SR(s).session = sessions(s);
        SR(s).runs = reshape(runs,length(runs),1);
    end
elseif ischar(runs)
    if strcmp(runs,'all')
        for s = 1:length(sessions)
            SR(s).session = sessions(s);
            if isempty(task)
                runs = str2double(spm_BIDS(BIDS,'runs','sub',subject,'ses',sprintf('%02.0f',sessions(s)),'type',scan_type,'modality',scan_modality));
            else
                runs = str2double(spm_BIDS(BIDS,'runs','sub',subject,'ses',sprintf('%02.0f',sessions(s)),'type',scan_type,'modality',scan_modality,'task',task));
            end
            SR(s).runs = reshape(runs,length(runs),1);
        end
    else
        error('In parameters: if ''runs'' is a char array it must be runs = ''all''')
    end
end

% Get filenames for every scan

for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    SR(s).filename = cell(length(runs),1);
    for r = 1:length(runs)
        run = runs(r);
        if isempty(task)
            filename = spm_BIDS(BIDS,'data','sub',subject,'ses',sprintf('%02.0f',session),'run',sprintf('%02.0f',run),'type',scan_type,'modality',scan_modality);
        else
            filename = spm_BIDS(BIDS,'data','sub',subject,'ses',sprintf('%02.0f',session),'run',sprintf('%02.0f',run),'task',task,'type',scan_type,'modality',scan_modality);
        end
        
        if strcmp(scan_type,'bold')
            if numel(filename) ~= 1; error('BIDS query did not select one & only one scan');end % this should never happen
            SR(s).filename{r} = filename{1};
            SR(s).namebase{r} = sprintf('sub-%s_ses-%02.0f_task-%s_run-%02.0f',subject,session,task,run);
        else
            SR(s).filename{r} = filename;
            [~,name,ext] = fileparts(filename{1});
            if strcmp(ext,'.gz'); [~,name] = fileparts(name); end
            stidx = strfind(name,['_' scan_type]);
            stl = length(['_' scan_type]);
            name(stidx:(stidx+stl-1)) = '';
            SR(s).namebase{r} = name;
        end
    end
end






