function [json_file,in_path] = get_corresp_json(input)

[in_path,in_name,ext] = fileparts(input);
if strcmp(ext,'.gz'); [in_path,in_name] = fileparts(in_name); end
json_file = fullfile(in_path,[in_name '.json']);