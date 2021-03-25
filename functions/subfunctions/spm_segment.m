function [tissue_files,tissues,seg_mat_files,mask_files] = spm_segment(T1_file,brain_prob,paths,T1_name,T2_file)
% script to use with human subjects
fprintf('Starting segmentation...\n')



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

%% T2
if ~exist('T2_file','var')
    T2 = '';
    zipT2 = 0;
else
    [pT2,T2name,ext] = fileparts(T2_file);
    zipT2 = 0;
    if strcmp(ext,'.gz')
        system(sprintf('gunzip %s',T2_file));
        T2_file = fullfile(pT2,T2name);
        [~,T2name] = fileparts(T2name);
        zipT2 = 1;
    end
    T2 = cellstr([T2_file ',1']);
end

%% Template
pathtemp = fullfile(fileparts(which('spm.m')),'tpm');
temp_file = fullfile(pathtemp,'TPM.nii');



%% matlabbatch
matlabbatch{1}.spm.spatial.preproc.channel(1).vols = T1;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel(1).write = [0 0];

if ~isempty(T2)
    matlabbatch{1}.spm.spatial.preproc.channel(2).vols = T2;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0];
end
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[temp_file ',1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[temp_file ',2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[temp_file ',3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[temp_file ',4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[temp_file ',5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[temp_file ',6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

%% Initialization
spm('defaults', 'FMRI');
spm_jobman('initcfg'); % initialization
spm_get_defaults('cmdline',true)

%% Run jobs
spm_jobman('run', matlabbatch)


%% Move files
if ~exist(paths.segmentation,'dir');mkdir(paths.segmentation);end % create segmentation folder if non-existant
movefile(fullfile(pT1,['c*' name '*']),paths.segmentation);
movefile(fullfile(pT1,[name '_seg8.mat']),paths.segmentation);


tissue_files = cell(3,1);
for i = 1:length(tissue_files)
    a = dir(fullfile(paths.segmentation,['c' num2str(i) '*' name '*']));
    tissue_files{i} = fullfile(paths.segmentation,a(1).name);
end
tissues = {'grey';'white';'csf'};


a = [name '_seg8.mat'];
if ~strcmp(name,T1_name)
    b = [T1_name a(length(name)+1:end)];
    movefile(fullfile(paths.segmentation,a),fullfile(paths.segmentation,b));
else
    b = a;
end
seg_mat_files{1} = fullfile(paths.segmentation,b);

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

if zipT2
    system(['gzip ' T2_file]);
end







%% Trash

% if ~exist('template','var') || isempty(template) || strcmp(template,'TPM')
%     pathtemp = fullfile(fileparts(which('spm.m')),'tpm');
%     temp_file = fullfile(pathtemp,'TPM.nii');
% elseif exist(template,'file')
%     [pTt,nTt,ext] = fileparts(template);
%     if strcmp(ext,'.gz')
%         system(sprintf('gunzip %s',template));
%         template = fullfile(pTt,nTt);
%         zipTt = 1;
%     else
%         zipTt = 0;
%     end
%     temp_file = template;
% else
%     error('Incorrect template')
% end