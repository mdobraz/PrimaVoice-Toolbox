function [tissue_files,tissues,seg_mat_files,mask_files] = spm_old_segment(T1_file,brain_prob,grey,white,csf,paths,T1_name)
% script to use with macaque subjects
% Note that multi-spectral (when there are two or more
% registered images of different contrasts) processing is
% not yet implemented for this method. (SPM manual)

fprintf('Starting SPM segmentation...\n')

%% T1
[pT1,name,ext] = fileparts(T1_file);
zipT1 = 0;
if strcmp(ext,'.gz')
    system(sprintf('gunzip %s',T1_file));
    T1_file = fullfile(pT1,name);
    [~,name] = fileparts(name);
    zipT1 = 1;
end
T1 = cellstr([T1_file ',1']);

%% Tissue Prob Maps
TPM = {grey;white;csf};
zipTPM = zeros(3,1);
for i = 1:length(TPM)
    [p,n,ext] = fileparts(TPM{i});
    if strcmp(ext,'.gz')
        system(sprintf('gunzip %s',TPM{i}));
        TPM{i} = fullfile(p,n);
        zipTPM(i) = 1;
    end
end

%% matlabbatch
matlabbatch{1}.spm.tools.oldseg.data = T1;
matlabbatch{1}.spm.tools.oldseg.output.GM = [0 0 1];
matlabbatch{1}.spm.tools.oldseg.output.WM = [0 0 1];
matlabbatch{1}.spm.tools.oldseg.output.CSF = [0 0 1];
matlabbatch{1}.spm.tools.oldseg.output.biascor = 0;
matlabbatch{1}.spm.tools.oldseg.output.cleanup = 0;
matlabbatch{1}.spm.tools.oldseg.opts.tpm = TPM;
matlabbatch{1}.spm.tools.oldseg.opts.ngaus = [2
                                              2
                                              2
                                              4];
matlabbatch{1}.spm.tools.oldseg.opts.regtype = ''; % 'subj'
matlabbatch{1}.spm.tools.oldseg.opts.warpreg = 1;
matlabbatch{1}.spm.tools.oldseg.opts.warpco = 25;
matlabbatch{1}.spm.tools.oldseg.opts.biasreg = 0.0001;
matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm = 60;
matlabbatch{1}.spm.tools.oldseg.opts.samp = 3;
matlabbatch{1}.spm.tools.oldseg.opts.msk = {''};

%% Initialization
spm('defaults', 'FMRI');
spm_jobman('initcfg'); % initialization
spm_get_defaults('cmdline',true)

%% Run jobs
spm_jobman('run', matlabbatch)



%% Move files
if ~exist(paths.segmentation,'dir');mkdir(paths.segmentation);end % create segmentation folder if non-existant
movefile(fullfile(pT1,['c*' name '*']),paths.segmentation);
movefile(fullfile(pT1,[name '_seg_*.mat']),paths.segmentation);


tissue_files = cell(3,1);
for i = 1:length(tissue_files)
    a = dir(fullfile(paths.segmentation,['c' num2str(i) '*' name '*']));
    tissue_files{i} = fullfile(paths.segmentation,a(1).name);
end
tissues = {'grey';'white';'csf'};


a = dir(fullfile(paths.segmentation,[name '_seg_*.mat']));
seg_mat_files = cell(length(a),1);
for i = 1:length(a)
    if ~strcmp(name,T1_name)
        b = [T1_name a(i).name(length(name)+1:end)];
        movefile(fullfile(paths.segmentation,a(i).name),fullfile(paths.segmentation,b));
    else
        b = a(i).name;
    end
    seg_mat_files{i} = fullfile(paths.segmentation,b);
end

%% Create brain masks from segmentation
mask_files = create_seg_masks(T1_file,brain_prob,tissue_files,tissues,paths);

if ~strcmp(name,T1_name)
    mffn = fieldnames(mask_files);
    for i = 1:length(mffn)
        [~,mname] = fileparts(fileparts2(mask_files.(mffn{i})));
        b = [T1_name mname(length(name)+1:end) '.nii.gz'];
        movefile(mask_files.(mffn{i}),fullfile(paths.segmentation,b));
        mask_files.(mffn{i}) = fullfile(paths.segmentation,b);
    end
end


%% Mask tissues
for i = 1:length(tissue_files)
    system(sprintf('%sfslmaths %s -mas %s %s',paths.FSL_prefix,tissue_files{i},mask_files.brain_mask,tissue_files{i}));
    delete(tissue_files{i});
    tissue_files{i} = [tissue_files{i} '.gz'];
end

%% Re-zip in-files if needed
if zipT1
    system(['gzip ' T1_file]);
end

for i = 1:length(TPM)
    if zipTPM(i)
        system(sprintf('gzip %s',TPM{i}));
    end
end



