function [xfm_i2r_file,xfm_r2i_file,flirt_out_file,brain_out_file,flirt_mask_out_file] = strong_register_FLIRT(in_file,flirt_in_file,ref_file,ref_brain_mask,cost,paths)

[in_path,in_name,ext] = fileparts(in_file);
if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end

[~,ref_name,ext] = fileparts(ref_file);
if strcmp(ext,'.gz'); [~,ref_name] = fileparts(ref_name); end


%% Calculate affine transform matrix from in_file to ref_file
reg_iterations = 3;

xfm_i2r_file = fullfile(in_path,sprintf('%s_FLIRT_to_%s.xfm',in_name,ref_name));
xfm_r2i_file = fullfile(in_path,sprintf('%s_FLIRT_to_%s.xfm',ref_name,in_name));
flirt_out_file = fullfile(in_path,sprintf('%s_FLIRT_to_%s.nii',in_name,ref_name));
if exist(flirt_out_file,'file') == 2; delete(flirt_out_file);end
flirt_mask_out_file = fullfile(in_path,[in_name '_brain_mask.nii']);
brain_out_file = fullfile(in_path,[in_name '_brain.nii']);
if exist(brain_out_file,'file') == 2; delete(brain_out_file);end
for i = 1:reg_iterations
    fprintf('FLIRT & mask iteration %i/%i...',i,reg_iterations)
    system(sprintf('%sflirt -in %s -ref %s -omat %s -out %s -dof 6 -cost %s -searchcost %s',paths.FSL_prefix,flirt_in_file,ref_file,xfm_i2r_file,flirt_out_file,cost,cost)); % compute xfm
	system(sprintf('%sconvert_xfm -omat %s -inverse %s',paths.FSL_prefix,xfm_r2i_file,xfm_i2r_file)); % inverse xfm
    system(sprintf('%sflirt -in %s -ref %s -out %s -interp nearestneighbour -applyxfm -init %s',paths.FSL_prefix,ref_brain_mask,flirt_in_file,flirt_mask_out_file,xfm_r2i_file)); % move brain mask to in_file
    system(sprintf('%sfslmaths %s -mul %s %s',paths.FSL_prefix,in_file,flirt_mask_out_file,brain_out_file)); % extract brain from mask
    
    flirt_in_file = brain_out_file;
    fprintf(' done.\n')
end

system(sprintf('gunzip %s',flirt_out_file));
system(sprintf('gunzip %s',brain_out_file));
[~,flirt_out_name] = fileparts(flirt_out_file);
view_slice_overlay(ref_file,flirt_out_file,0)
figurewrite(fullfile(in_path,flirt_out_name),[],[],[],1) % Using GLMdenoise function

return

