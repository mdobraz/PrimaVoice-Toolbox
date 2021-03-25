function [json,Pstep] = get_preproc_history(input)

json_file = get_corresp_json(input);

if exist(json_file,'file')
    json = spm_jsonread(json_file);
    fn = fieldnames(json);
    IndexC = strfind(fn,'PV_preprocessing_step');
    Index = not(cellfun('isempty',IndexC));
    if any(Index)
        Ps = sort(fn(Index));
        Pstep = str2double(Ps{end}(end-1:end));
    else
        json = [];
        Pstep = 0;
    end
else
    json = [];
    Pstep = 0;
end



