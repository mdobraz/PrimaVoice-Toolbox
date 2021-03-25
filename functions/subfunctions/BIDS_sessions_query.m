function result = BIDS_sessions_query(BIDS,subject,modality,type,task)

nr = 0;
for s = 1:numel(BIDS.subjects)
    if strcmp(['sub-' subject],BIDS.subjects(s).name)
        scans = BIDS.subjects(s).(modality);
        for scan = 1:numel(scans)
            res = BIDS.subjects(s).(modality)(scan).ses;
            if exist('type','var')
                if ~strcmp(type,BIDS.subjects(s).(modality)(scan).type)
                    res = '';
                end
            end
            
            if exist('task','var')
                if ~strcmp(task,BIDS.subjects(s).(modality)(scan).task)
                    res = '';
                end
            end
            if ~isempty(res)
                nr = nr + 1;
                result{nr} = res;
            end
        end
    end
end

result = unique(result);