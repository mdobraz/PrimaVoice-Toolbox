function [base,ext]=fileparts2(string)
[path,name,ext]=fileparts(deblank(string));
if strcmp(ext,'.gz')
   [~,name,ext] = fileparts(name);
end   
if isempty(path)
   base = name;
else
   base = fullfile(path,name);
end
return

