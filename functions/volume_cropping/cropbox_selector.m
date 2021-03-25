%% Choose file to load
startdir = fileparts(file_to_load);

[CB_to_load,path_cbtl] = uigetfile(fullfile(startdir,'*.cropbox'),'Select file to load');

if CB_to_load == 0
    error('Please select one file to load')
end

CB_to_load = [path_cbtl CB_to_load];
