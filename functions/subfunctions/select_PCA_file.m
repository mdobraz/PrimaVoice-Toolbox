function [mat_file,nii_file] = select_PCA_file(infile,paths,mask_name)

if ~exist('mask_name','var') || isempty(mask_name)
    mask_name = 'automask';
end

[~,name] = fileparts(infile);
out_base = fullfile(paths.PCA,[name '_PCA_' mask_name]);

mat_file = [out_base '.mat'];
nii_file = [out_base '_pca.nii'];
