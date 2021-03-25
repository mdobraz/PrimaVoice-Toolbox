function out_file = spm_sanlm(in_file,out_dir)

fprintf('\nStarting denoising...\n\n')

[in_path,in_name,ext] = fileparts(in_file);

if strcmp(ext,'.gz')
    [~,in_name,ext] = fileparts(in_name);
	in_file = fullfile(in_path,[in_name ext]);
	if exist(in_file,'file') == 2; delete(in_file);end
	system(sprintf('gunzip %s',[in_file '.gz']));
    tozip = 1;
else
    tozip = 0;
end

spm('defaults', 'FMRI');
spm_jobman('initcfg'); % initialization
spm_get_defaults('cmdline',true)

spm_prefix = 'sanlm_';
out_file = fullfile(in_path,[in_name '_denoised' ext]);

matlabbatch{1}.spm.tools.cat.tools.sanlm.prefix = spm_prefix;
matlabbatch{1}.spm.tools.cat.tools.sanlm.NCstr = Inf;
matlabbatch{1}.spm.tools.cat.tools.sanlm.rician = 0;


spm_out_file = fullfile(in_path,[spm_prefix in_name ext]);

matlabbatch{1}.spm.tools.cat.tools.sanlm.data = {in_file};
spm_jobman('run', matlabbatch);
movefile(spm_out_file,out_file);

if tozip
    system(['gzip ' in_file]);
end

system(['gzip ' out_file]);
out_file = [out_file '.gz'];


if exist('out_dir','var') % move file if a out_dir has been specified
    [~,out_name,ext] = fileparts(out_file);
    out_nopath = [out_name ext];
    out_file_in_dir = fullfile(out_dir,out_nopath);
    movefile(out_file,out_file_in_dir);
    out_file = out_file_in_dir;
end




