function apply_extent_threshold(in_file,out_file,clust_thres,extent_threshold_vox,isMION,conn)

if ~exist('conn','var')
	conn = 18;
end

if ~exist('isMION','var')
	isMION = 0;
end

[~,~,ext] = fileparts(out_file);
if strcmp(ext,'.gz')
	out_file = out_file(1:end-3);
    tozip = 1;
else
	tozip = 0;
end


Pt = spm_vol(in_file);
Yt = spm_read_vols(Pt);

if isMION
    Yt = -Yt;
end


%% extent threshold
% mask using the cluster threshold
mask = false(size(Yt));
mask(Yt > clust_thres) = true;

% remove clusters smaller than extent threshold in mask
[labels,num] = bwlabeln(mask,conn); % 6-connected, 18-connected or 26-connected
for n = 1:num
    cluster_size = sum(labels(:) == n);
    if cluster_size < extent_threshold_vox
        mask(labels == n) = false;
    end
end

Yt(~mask) = 0;

P = Pt;
P.fname = out_file;
P.descrip = 'Cluster thresholded';
spm_write_vol(P,Yt);

if tozip
	system(['gzip ' out_file]);
end


% in_file = '/hpc/banco/Primavoice_Data_and_Analysis/second_level_analysis_macaques/second_level_macaques_fixed_macaque_vs_all_t.nii.gz';
% out_file = '/hpc/banco/Primavoice_Data_and_Analysis/second_level_analysis_macaques/second_level_macaques_fixed_macaque_vs_all_t_clust-thresholded.nii';