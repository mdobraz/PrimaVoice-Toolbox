%% Settings, path...
if exist(paths.segmentation,'dir')
    system(['rm -r ' paths.segmentation]);
else
    mkdir(paths.segmentation); % create segmentation folder if non-existant
end
if ~exist(paths.preprocessing,'dir');mkdir(paths.preprocessing);end % create folder if non-existant
if ~exist(paths.seg_fig,'dir');mkdir(paths.seg_fig);end % create folder if non-existant

setenv('ANTSPATH',paths.ANTS_path);

T2_exists = 0;
if isfield(paths,'T2_anat_suffix') && ~isempty(paths.T2_anat_suffix)
    T2_exists = 1;
end

% Setting paths & -p option for bash scripts
debias_path = which('T1xT2BiasFieldCorrection.sh');
BET_path = which('T1xT2BET.sh');
IterREGBET_path = which('IterREGBET.sh');
CropVolume_path = which('CropVolume.sh');
if isfield(paths,'FSL_prefix') && ~isempty(paths.FSL_prefix)
    popt = ['-p ' paths.FSL_prefix];
else
    popt ='';
end



%% Denoising
T1_denoised = spm_sanlm(T1w_files.raw,paths.preprocessing);
if T2_exists
    T2_denoised = spm_sanlm(T2w_files.raw,paths.preprocessing);
    aT2opt = '-aT2';
else
    T2_denoised = T1_denoised;
    aT2opt = '';
end



%% First large Debias (brain_seg.BETf * 0.6)
system(sprintf('bash %s -t1 %s -t2 %s %s -bet 2 -f %s -g %s %s',debias_path,T1_denoised,T2_denoised,aT2opt,num2str(brain_seg.BETf * 0.6),num2str(brain_seg.BETg),popt));
T1_debiased = [fileparts2(T1_denoised) '_debiased.nii.gz'];
T2_debiased = [fileparts2(T2_denoised) '_debiased.nii.gz'];
if T2_exists
    T2_denoised_in_T1 =  [fileparts2(T2_denoised) '-in-T1w.nii.gz'];
    movefile(T2_denoised_in_T1,T2_denoised)
end
T1_BET_mask = [fileparts2(T1_debiased) '_BET_mask.nii.gz'];
view_slice_overlay(T1_debiased,T1_BET_mask,0,[],[],[],lines(2))
figurewrite(fileparts2(T1_BET_mask),[],[],[],1) % Using GLMdenoise function

%% T1xT2BET
system(sprintf('bash %s -t1 %s -t2 %s -n 3 -f %s -g %s %s',BET_path,T1_debiased,T2_debiased,num2str(brain_seg.BETf),num2str(brain_seg.BETg),popt));
T1_BET = [fileparts2(T1_debiased) '_BET.nii.gz'];
% save figure
rel_cog = vol_rel_cog(T1_BET,paths);
view_slice_overlay(T1_debiased,T1_BET,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(T1_BET),[],[],[],1) % Using GLMdenoise function

%% IterREGBET for non-humans
if ~strcmp(paths.species,'human')
    system(sprintf('bash %s -inw %s -inb %s -refb %s -cost %s -dof 12 -n 3 %s',IterREGBET_path,T1_debiased,T1_BET,paths.template_brain,brain_seg.flirt_cost,popt));
    T1_BET = [fileparts2(T1_debiased) '_IRbrain.nii.gz'];
    % save figure
    rel_cog = vol_rel_cog(T1_BET,paths);
    view_slice_overlay(T1_debiased,T1_BET,0,[],[],[],[],rel_cog)
    figurewrite(fileparts2(T1_BET),[],[],[],1) % Using GLMdenoise function
end




if brain_seg.redebias % crop & re-run everything with previously extracted brain
    %% Crop files
    system(sprintf('bash %s -i %s -i %s -i %s -b %s -c 20 %s',CropVolume_path,T1_denoised,T2_denoised,T1_BET,T1_BET,popt));
    T1_cropped = [fileparts2(T1_denoised) '_cropped.nii.gz'];
    T2_cropped = [fileparts2(T2_denoised) '_cropped.nii.gz'];
    T1_BET_cropped = [fileparts2(T1_BET) '_cropped.nii.gz'];
    T1_cropped_BET = [fileparts2(T1_denoised) '_cropped_BET.nii.gz'];
    movefile(T1_BET_cropped,T1_cropped_BET)
    %%%% Rename file


    %% Debias with brain extracted
    system(sprintf('bash %s -t1 %s -t2 %s -b %s %s',debias_path,T1_cropped,T2_cropped,T1_cropped_BET,popt));
    T1_debiased = [fileparts2(T1_cropped) '_debiased.nii.gz'];
    T2_debiased = [fileparts2(T2_cropped) '_debiased.nii.gz'];


    %% T1xT2BET
    system(sprintf('bash %s -t1 %s -t2 %s -n 3 -f %s %s',BET_path,T1_debiased,T2_debiased,num2str(brain_seg.BETf),popt));
    T1_BET = [fileparts2(T1_debiased) '_BET.nii.gz'];
    % save figure
    rel_cog = vol_rel_cog(T1_BET,paths);
    view_slice_overlay(T1_debiased,T1_BET,0,[],[],[],[],rel_cog)
    figurewrite(fileparts2(T1_BET),[],[],[],1) % Using GLMdenoise function


    %% IterREGBET for non-humans
    if ~strcmp(paths.species,'human')
        system(sprintf('bash %s -inw %s -inb %s -refb %s -cost %s -dof 12 -n 3 %s',IterREGBET_path,T1_debiased,T1_BET,paths.template_brain,brain_seg.flirt_cost,popt));
        T1_BET = [fileparts2(T1_debiased) '_IRbrain.nii.gz'];
        % save figure
        rel_cog = vol_rel_cog(T1_BET,paths);
        view_slice_overlay(T1_debiased,T1_BET,0,[],[],[],[],rel_cog)
        figurewrite(fileparts2(T1_BET),[],[],[],1) % Using GLMdenoise function
    end
end



%% Rigid body registration in order to put the T1 & T2 in the template space
T1_brain_in_temp = [fileparts2(T1_BET) '_in-' paths.template_name '.nii.gz'];
anat2temp_xfm = fullfile(paths.preprocessing,'T1_rigid_to_template.xfm');
system(sprintf('%sflirt -in %s -ref %s -dof 6 -cost %s -out %s -omat %s',paths.FSL_prefix,T1_BET,paths.template_brain,brain_seg.flirt_cost,T1_brain_in_temp,anat2temp_xfm));

[~,T1_debiased_name,ext] = fileparts(T1_debiased);
if strcmp(ext,'.gz'); [~,T1_debiased_name,ext] = fileparts(T1_debiased_name); end
T1_in_temp = fullfile(paths.anat_dir,[T1_debiased_name '_in-' paths.template_name '.nii.gz']);
system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s',paths.FSL_prefix,T1_debiased,paths.template_brain,T1_in_temp,anat2temp_xfm));

[~,T2_debiased_name,ext] = fileparts(T2_debiased);
if strcmp(ext,'.gz'); [~,T2_debiased_name,ext] = fileparts(T2_debiased_name); end
T2_in_temp = fullfile(paths.anat_dir,[T2_debiased_name '_in-' paths.template_name '.nii.gz']);
system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s',paths.FSL_prefix,T2_debiased,paths.template_brain,T2_in_temp,anat2temp_xfm));

T2_brain_in_temp = [fileparts2(T2_debiased) '_BET_in-' paths.template_name '.nii.gz'];
system(sprintf('%sfslmaths %s -mas %s %s -odt short',paths.FSL_prefix,T2_in_temp,T1_brain_in_temp,T2_brain_in_temp));



%% ANTs Registration (affine + non-linear)
if brain_seg.longSyN
    aRS = fullfile(paths.ANTS_scripts,'antsRegistrationSyN.sh');
else
    aRS = fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh');    
end
system(sprintf('%s -d 3 -f %s -f %s -m %s -m %s -o %s -j 1',...
    aRS,T1_brain_in_temp,T1_in_temp,...
    paths.template_brain,paths.template_head,paths.ANTs_out_base));

paths.temp_to_anat_xfm = [paths.ANTs_out_base '0GenericAffine.mat'];
paths.anat_to_temp_xfm = [paths.ANTs_out_base '0InverseGenericAffine.mat'];
paths.temp_to_anat_warp = [paths.ANTs_out_base '1Warp.nii.gz'];
paths.anat_to_temp_warp = [paths.ANTs_out_base '1InverseWarp.nii.gz'];

%% Invert affine
c3d_path = which('c3d_affine_tool.b');
system(sprintf('%s -itk %s -inv -oitk %s',c3d_path,paths.temp_to_anat_xfm,paths.anat_to_temp_xfm));

%% Apply transforms to tissues
% the order of transformations has to be inverted in the antsApplyTransforms call

fprintf('Applying warps...\n')
to_move = {'paths.template_tissue.grey';'paths.template_tissue.white';'paths.template_tissue.csf'};
reg_prob_maps = cell(length(to_move),1);
for i = 1:length(to_move)
    in_file = eval(to_move{i});
    [~,in_name,ext] = fileparts(in_file);
    if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end
    reg_prob_maps{i} = fullfile(paths.segmentation,[in_name '-in-anat.nii.gz']);
    system(sprintf('%s -i %s -r %s -o %s -t %s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),in_file,T1_in_temp,reg_prob_maps{i},paths.temp_to_anat_warp,paths.temp_to_anat_xfm));
end


%% Choose image to segment (brain extracted or not)
if brain_seg.use_brain_for_seg
    T1_to_seg = T1_brain_in_temp;
    if T2_exists
        T2_to_seg = T2_brain_in_temp;
    end
else
    T1_to_seg = T1_in_temp;
    if T2_exists
        T2_to_seg = T2_in_temp;
    end
end
[~,T1_in_temp_name] = fileparts(fileparts2(T1_in_temp));


%% Segment
if strcmp(paths.species,'human') && brain_seg.use_new_segment
    if T2_exists
        [tissue_files,tissues,seg_mat_files,mask_files] = spm_segment(T1_to_seg,brain_seg.brain_mask_prob,paths,T1_in_temp_name,T2_to_seg);
    else
        [tissue_files,tissues,seg_mat_files,mask_files] = spm_segment(T1_to_seg,brain_seg.brain_mask_prob,paths,T1_in_temp_name);
    end
else
    [tissue_files,tissues,seg_mat_files,mask_files] = spm_old_segment(T1_to_seg,brain_seg.brain_mask_prob,reg_prob_maps{1},reg_prob_maps{2},reg_prob_maps{3},paths,T1_in_temp_name);
end



%% Re-register if wanted & perform segmentation again
if brain_seg.rereg
    system(sprintf('%sfslmaths %s -mas %s %s',paths.FSL_prefix,T1_in_temp,mask_files.brain_mask,T1_brain_in_temp));
    for i = 1:length(tissue_files)
        delete(tissue_files{i});
    end
    mffn = fieldnames(mask_files);
    for i = 1:length(mffn)
        delete(mask_files.(mffn{i}));
    end

    fprintf('Applying warps...\n')
    to_move = {'paths.template_brain';'paths.template_head'};
    reg_temp = cell(length(to_move),1);
    for i = 1:length(to_move)
        in_file = eval(to_move{i});
        [~,in_name,ext] = fileparts(in_file);
        if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end
        reg_temp{i} = fullfile(paths.segmentation,[in_name '-in-anat.nii.gz']);
        system(sprintf('%s -i %s -r %s -o %s -t %s -t %s -n Linear',fullfile('$ANTSPATH','antsApplyTransforms'),in_file,T1_in_temp,reg_temp{i},paths.temp_to_anat_warp,paths.temp_to_anat_xfm));
    end

    ANTSout = fullfile(paths.segmentation,'reSyN_template_to_anat');
    system(sprintf('%s -d 3 -f %s -f %s -m %s -m %s -o %s -j 1',...
        aRS,T1_brain_in_temp,T1_in_temp,...
        reg_temp{1},reg_temp{2},ANTSout));

    re_temp_to_anat_xfm = [ANTSout '0GenericAffine.mat'];
    re_temp_to_anat_warp = [ANTSout '1Warp.nii.gz'];

    fprintf('Applying warps...\n')
    re_reg_prob_maps = cell(length(reg_prob_maps),1);
    for i = 1:length(reg_prob_maps)
        in_file = reg_prob_maps{i};
        [~,in_name,ext] = fileparts(in_file);
        if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end
        re_reg_prob_maps{i} = fullfile(paths.segmentation,['re_' in_name '.nii.gz']);
        system(sprintf('%s -i %s -r %s -o %s -t %s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),in_file,T1_in_temp,re_reg_prob_maps{i},re_temp_to_anat_warp,re_temp_to_anat_xfm));
    end

    if strcmp(paths.species,'human') && brain_seg.use_new_segment
        if T2_exists
            [tissue_files,tissues,seg_mat_files,mask_files] = spm_segment(T1_to_seg,brain_seg.brain_mask_prob,paths,T1_name,T2_to_seg);
        else
            [tissue_files,tissues,seg_mat_files,mask_files] = spm_segment(T1_to_seg,brain_seg.brain_mask_prob,paths,T1_in_temp_name);
        end
    else
        [tissue_files,tissues,seg_mat_files,mask_files] = spm_old_segment(T1_to_seg,brain_seg.brain_mask_prob,reg_prob_maps{1},reg_prob_maps{2},reg_prob_maps{3},paths,T1_in_temp_name);
    end
end

%% save figures
warped_temp = [paths.ANTs_out_base 'Warped.nii.gz'];
rel_cog = vol_rel_cog(warped_temp,paths);
view_slice_overlay(T1_in_temp,warped_temp,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(warped_temp),[],[],[],1) % Using GLMdenoise function

warped_anat = [paths.ANTs_out_base 'InverseWarped.nii.gz'];
rel_cog = vol_rel_cog(warped_anat,paths);
view_slice_overlay(paths.template_head,warped_anat,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(warped_anat),[],[],[],1) % Using GLMdenoise function

rel_cog = vol_rel_cog(mask_files.brain,paths);
rel_cog(1) = rel_cog(1)*1.5;
view_slice_overlay(T1_in_temp,mask_files.brain_segmented,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(mask_files.brain_segmented),[],[],[],1) % Using GLMdenoise function



%% T1w_files
T1w_files.(paths.template_name).name = T1_in_temp_name;
T1w_files.(paths.template_name).fullhead = T1_in_temp;
T1w_files.(paths.template_name).brain = mask_files.brain;
T1w_files.(paths.template_name).brain_mask = mask_files.brain_mask;
T1w_files.(paths.template_name).brain_segmented = mask_files.brain_segmented;
T1w_files.(paths.template_name).brain_seg_mat_files = seg_mat_files;
T1w_files.(paths.template_name).segmentation_parameters = brain_seg;

if T2_exists
    [~,T2_in_temp_name] = fileparts(fileparts2(T2_in_temp));
    T2w_files.(paths.template_name).name = T2_in_temp_name;
    T2w_files.(paths.template_name).fullhead = T2_in_temp;
end


fprintf('Creating prob masks...')
%% Create 01, 50 & 99 percents masks for each tissue
probs = paths.anat.tissue_probs;
for i = 1:numel(tissues)
    T1w_files.(paths.template_name).tissue_probs.(tissues{i}) = fullfile(paths.segmentation,sprintf('%s_%s_prob.nii.gz',T1_in_temp_name,tissues{i}));
    
    for p = 1:numel(probs)
        % Create a mask of voxels with prob >= probs-percent
        T1w_files.(paths.template_name).tissue_bins.(tissues{i}).(sprintf('p%s',probs{p})) = fullfile(paths.segmentation,sprintf('%s_%s_p%s.nii.gz',T1_in_temp_name,tissues{i},probs{p}));
        maths_out_file = T1w_files.(paths.template_name).tissue_bins.(tissues{i}).(sprintf('p%s',probs{p}));
        system(sprintf('%sfslmaths %s -thr 0.%s %s',paths.FSL_prefix,tissue_files{i},probs{p},maths_out_file));
        system(sprintf('%sfslmaths %s -bin %s -odt short',paths.FSL_prefix,maths_out_file,maths_out_file));
        
        % Keep the largest cluster
        system(['gunzip ' maths_out_file]);
        maths_out_file(end-2:end) = '';
        P = spm_vol(maths_out_file);
        Y = spm_read_vols(P);
        Y = keep_largest_cluster(Y,6);
        spm_write_vol(P,Y);
        system(['gzip ' maths_out_file]);
        maths_out_file = [maths_out_file '.gz'];
        % system(sprintf('%s 3 %s GetLargestComponent %s',fullfile('$ANTSPATH','ImageMath'),maths_out_file,maths_out_file));
        system(sprintf('%sfslmaths %s %s -odt short',paths.FSL_prefix,maths_out_file,maths_out_file)); % convert to short

        %% Create surface
        % if ~exist(paths.surfaces,'dir');mkdir(paths.surfaces);end % create folder if non-existant
        %     [~,mof_name] = fileparts(fileparts2(maths_out_file));
        % system(sprintf('%s -i %s',fullfile(paths.brainvisa,'AimsMesh'),fullfile(paths.segmentation,mof_name))); % Mesh
        % gii = [mof_name '_1_0.gii'];
        % movefile(fullfile(paths.segmentation,gii),paths.surfaces);
        % system(sprintf('%s -i %s --algoType laplacian -o %s',fullfile(paths.brainvisa,'AimsMeshSmoothing'),fullfile(paths.surfaces,gii),fullfile(paths.surfaces,['smoothed_' gii]))); % Smooth
    end
    movefile(tissue_files{i},T1w_files.(paths.template_name).tissue_probs.(tissues{i}));
end


%% Create white + csf masks for PCA
pca_tissues = {'white';'csf'};
grey_mask = T1w_files.(paths.template_name).tissue_bins.(tissues{1}).(sprintf('p%s',probs{1}));
for i = 1:numel(pca_tissues)
    in_file = T1w_files.(paths.template_name).tissue_bins.(pca_tissues{i}).(sprintf('p%s',probs{end}));
    out_file = fullfile(paths.segmentation,sprintf('%s_%s_pca_mask.nii.gz',T1_in_temp_name,pca_tissues{i}));

    system(sprintf('%sfslmaths %s -sub %s %s',paths.FSL_prefix,in_file,grey_mask,out_file)); % remove voxels that might be in the grey matter
    system(sprintf('%sfslmaths %s -bin %s',paths.FSL_prefix,out_file,out_file)); % binarize
    system(sprintf('%sfslmaths %s -kernel boxv 3 -ero %s',paths.FSL_prefix,out_file,out_file)); % erode mask (to be sure no voxels of interest will be selected)

    T1w_files.(paths.template_name).pca_masks.(pca_tissues{i}) = out_file;
    % save figure
    view_slice_overlay(T1_in_temp,out_file,0,[],[],[],lines(2),rel_cog)
    figurewrite(fileparts2(out_file),[],[],[],1) % Using GLMdenoise function
end




movefile(fullfile(paths.preprocessing,'*.png'),paths.seg_fig);
movefile(fullfile(paths.segmentation,'*.png'),paths.seg_fig);


%% Save T1w_files & T2w_files
spm_jsonwrite(T1w_files_json,T1w_files,struct('indent','    '));
if T2_exists
    spm_jsonwrite(T2w_files_json,T2w_files,struct('indent','    '));
end

fprintf(' done.\n')


%% Trash

%% Create brain mask from concatenation of 3 tissues at 1% prob
% T1_brainmask_in_temp = fullfile(paths.segmentation,[T1w_files.name '-brainmask-in-template.nii.gz']);
% system(sprintf('%sfslmaths %s -add %s -add %s -bin %s -odt short',paths.FSL_prefix,fullfile(paths.segmentation,[T1w_files.name '_white_01p.nii.gz']),fullfile(paths.segmentation,[T1w_files.name '_grey_01p.nii.gz']),fullfile(paths.segmentation,[T1w_files.name '_csf_01p.nii.gz']),T1_brainmask_in_temp));
% system(sprintf('%s 3 %s FillHoles %s',fullfile('$ANTSPATH','ImageMath'),T1_brainmask_in_temp,T1_brainmask_in_temp));
% system(sprintf('%sfslmaths %s %s -odt short',paths.FSL_prefix,T1_brainmask_in_temp,T1_brainmask_in_temp)); % convert to short

% % mask brain
% system(sprintf('%sfslmaths %s -mas %s %s -odt short',paths.FSL_prefix,T1_in_temp,T1_brainmask_in_temp,T1_brain_in_temp));


