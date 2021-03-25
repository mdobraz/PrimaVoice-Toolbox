function [PCs,EVs] = run_PCA_image(infile,paths,mask_name,mask_file)

if ~exist('mask_name','var') || isempty(mask_name)
    mask_name = 'automask';
end

if ~exist('mask_file','var') || isempty(mask_file)
    mask_file = infile;
end

if ~exist(paths.PCA,'dir');mkdir(paths.PCA);end % create folder if non-existant
nPCs = 20;

[~,name] = fileparts(infile);
out_base = fullfile(paths.PCA,[name '_PCA_' mask_name]);
figureprep([0 0 1600 800]);
[PCs,~,EVs] = pca_image(infile,[],nPCs,mask_file,[],out_base);
figurewrite(out_base,[],0); % the 0 is to force eps figure
save(out_base,'PCs','EVs')

