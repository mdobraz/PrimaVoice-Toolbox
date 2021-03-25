function coincidence_map(t_files,paths,peak_thres,ext_thres,clust_thres,out_base,inverse)

fprintf('\nComputing coincidence maps...')

IDs = (1:numel(t_files)).^2;

for i = 1:numel(t_files)
    % load t-file
    Pt = spm_vol(t_files{i});
    Yt = spm_read_vols(Pt);

    % if MION, inverse values
    if strcmpi(paths.contrast_agent,'MION')
        Yt = -Yt;
    end

    % if inverse, inverse values
    if exist('inverse','var') && inverse
        Yt = -Yt;
    end

    % create mask
    if ~exist('maskP','var')
        maskP = zeros(size(Yt));
        maskC = zeros(size(Yt));
        maskPa = zeros(size(Yt));
        maskCa = zeros(size(Yt));
    end

    %% peak threshold
    maskP(Yt >= peak_thres(i)) = maskP(Yt >= peak_thres(i)) + 1;
    maskPa(Yt >= peak_thres(i)) = maskPa(Yt >= peak_thres(i)) + IDs(i);

    %% extent threshold
    % mask using the cluster threshold
    mask = false(size(Yt));
    mask(Yt > clust_thres(i)) = true;

    hd = fmris_read_image(t_files{i});
    extent_threshold_vox = round(ext_thres(i) / abs(prod(hd.vox)));

    % remove clusters smaller than extent threshold in mask
    [labels,num] = bwlabeln(mask,18); % 6-connected, 18-connected or 26-connected
    for n = 1:num
        cluster_size = sum(labels(:) == n);
        if cluster_size < extent_threshold_vox
            mask(labels == n) = false;
        end
    end
    maskC(mask) = maskC(mask) + 1;
    maskCa(mask) = maskCa(mask) + IDs(i);

    
end

P = Pt;
P.fname = [out_base '_peak_coincidence_map.nii'];
P.descrip = 'Peak occurences';
spm_write_vol(P,maskP);

P.fname = [out_base '_cluster_coincidence_map.nii'];
P.descrip = 'Cluster occurences';
spm_write_vol(P,maskC);

P.fname = [out_base '_peak_coincidence_map_detailled.nii'];
P.descrip = 'Peak occurences';
spm_write_vol(P,maskPa);

P.fname = [out_base '_cluster_coincidence_map_detailled.nii'];
P.descrip = 'Cluster occurences';
spm_write_vol(P,maskCa);

fprintf(' done.\n\n')
