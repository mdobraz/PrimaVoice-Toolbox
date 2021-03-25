function mask_files = create_seg_masks(T1_file,prob,tissue_files,tissues,paths)
% Create brain masks from segmentation files
fprintf('Creating masks...')



%% Load segmentation prob maps
for i = 1:numel(tissues)
    P.(tissues{i}) = spm_vol(tissue_files{i});
    Y.(tissues{i}) = spm_read_vols(P.(tissues{i}));
    Y.(tissues{i}) = keep_largest_cluster(Y.(tissues{i}),6,prob); %% keep largest clusters > prob, 6-connected
end


%% Create brain mask
brain_mask = Y.(tissues{1}) > prob | Y.(tissues{2}) > prob | Y.(tissues{3}) > prob; % Union of all tissues


%% Dilate & smooth mask
% dilated_brain_mask = erode_or_dilate(brain_mask,'dilate',18);
% smoothed_brain_mask = erode_or_dilate(dilated_brain_mask,'erode',18);

%% Fill holes
for i = 1:size(brain_mask,1)
    brain_mask(i,:,:) = imfill(brain_mask(i,:,:),'holes');
end
for j = 1:size(brain_mask,2)
    brain_mask(:,j,:) = imfill(brain_mask(:,j,:),'holes');
end
for k = 1:size(brain_mask,3)
    brain_mask(:,:,k) = imfill(brain_mask(:,:,k),'holes');
end

%% Keep largest cluster, 6-connected
brain_mask = keep_largest_cluster(brain_mask,6);


%% Map of the 3 tissues
Yall = zeros(size(Y.(tissues{1}),1),size(Y.(tissues{1}),2),size(Y.(tissues{1}),3),numel(tissues));
for i = 1:numel(tissues)
    Yall(:,:,:,i) = Y.(tissues{i});
end

Ymerged = zeros(size(Y.(tissues{1})));
for i = 1:size(Y.(tissues{1}),1)
    for j = 1:size(Y.(tissues{1}),2)
        for k = 1:size(Y.(tissues{1}),3)
            if brain_mask(i,j,k)
                [maxtissue,itissue] = max(Yall(i,j,k,:));
                Ymerged(i,j,k) = itissue;
            end
        end
    end
end


%% Keep only largest clusters of grey & white in merged, to affine brain mask
brain_mask = logical(Ymerged);
brain_mask(Ymerged~=3) = 0; % only keep csf
for i = 1:2
    X = Ymerged;
    X(Ymerged~=i) = 0;
    X = keep_largest_cluster(X,6);
    brain_mask((X~=0)) = 1;
end

%% Fill holes
for i = 1:size(brain_mask,1)
    brain_mask(i,:,:) = imfill(brain_mask(i,:,:),'holes');
end
for j = 1:size(brain_mask,2)
    brain_mask(:,j,:) = imfill(brain_mask(:,j,:),'holes');
end
for k = 1:size(brain_mask,3)
    brain_mask(:,:,k) = imfill(brain_mask(:,:,k),'holes');
end

Ymerged(~brain_mask) = 0;


%% save files
[~,name,ext] = fileparts(T1_file);

mask_files.brain = fullfile(paths.segmentation,[name '_brain.nii.gz']);
mask_files.brain_mask = fullfile(paths.segmentation,[name '_brain_mask' ext]);
mask_files.brain_segmented = fullfile(paths.segmentation,[name '_brain_segmented' ext]);

P = spm_vol(T1_file);

P.fname = mask_files.brain_mask;
spm_write_vol(P,brain_mask);

P.fname = mask_files.brain_segmented;
spm_write_vol(P,Ymerged);


% smooth brain mask
system(sprintf('%sfslmaths %s -kernel boxv 3 -ero %s -odt short',paths.FSL_prefix,mask_files.brain_mask,mask_files.brain_mask)); % erode mask 
delete(mask_files.brain_mask);
mask_files.brain_mask = [mask_files.brain_mask '.gz'];
system(sprintf('%sfslmaths %s -kernel boxv 3 -dilM %s -odt short',paths.FSL_prefix,mask_files.brain_mask,mask_files.brain_mask)); % dilate mask 


% mask sementation
system(sprintf('%sfslmaths %s -mas %s %s -odt short',paths.FSL_prefix,mask_files.brain_segmented,mask_files.brain_mask,mask_files.brain_segmented));
delete(mask_files.brain_segmented);
mask_files.brain_segmented = [mask_files.brain_segmented '.gz'];

% mask brain
system(sprintf('%sfslmaths %s -mas %s %s -odt short',paths.FSL_prefix,T1_file,mask_files.brain_mask,mask_files.brain));





fprintf(' done.\n')



%% Trash

% P.fname = fullfile(paths.segmentation,[name '_brain_mask' ext]);
% spm_write_vol(P,smoothed_brain_mask);

% P.fname = fullfile(paths.segmentation,[name '_brain_mask_dil' ext]);
% spm_write_vol(P,dilated_brain_mask);