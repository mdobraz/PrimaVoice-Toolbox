function file_to_load = file_selector(file_type,paths)
% Choose file to load

if ~exist('file_type','var') || isempty(file_type)
    file_type = '';
end

startdir = '';
if exist('paths','var')
    if isfield(paths,'anat_file')
        startdir = fileparts(paths.anat_file);
    else
        startdir = fileparts(paths);
    end
end
[file_to_load,path_ftl] = uigetfile(fullfile(startdir,['*' file_type '*']),'Select file to load');

if file_to_load == 0
    error('Please select one file to load')
end

file_to_load = [path_ftl file_to_load];
