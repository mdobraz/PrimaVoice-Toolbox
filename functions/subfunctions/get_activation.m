function [n_vox_mask,n_vox_maskP,meanT,meanTP,maxT] = get_activation(t_map_file,con_name,desired_p_val,p_val_peak,peak_threshold,extent_threshold_vox,Cluster_threshold,paths,saveit)

peak_threshold = peak_threshold(p_val_peak == desired_p_val);

%% load t-map
P = spm_vol(t_map_file);
Y = spm_read_vols(P);

%% if MION, inverse t-map
if strcmpi(paths.contrast_agent,'MION')
    Y = -Y;
end

%% create masks
maskP = false(size(Y));
maskC = false(size(Y));

%% mask using the peak threshold
maskP(Y > peak_threshold) = true;

%% mask using the cluster threshold
maskCtmp = maskC;
maskCtmp(Y > Cluster_threshold) = true;

% remove clusters smaller than extent threshold in maskCtmp
[labels,num] = bwlabeln(maskCtmp,18); % 6-connected, 18-connected or 26-connected
for n = 1:num
    cluster_size = sum(labels(:) == n);
    if cluster_size >= extent_threshold_vox
        maskC(labels == n) = true;
    end
end

mask = maskP;
mask(maskC) = true;


%% get activation in mask
n_vox_mask = sum(mask(:));
n_vox_maskP = sum(maskP(:));
meanT = mean(Y(mask>0));
meanTP = mean(Y(maskP>0));
maxT = max(Y(:));


%% save
if ~exist('saveit','var') || saveit
	if ~exist(paths.AC_folder,'dir');mkdir(paths.AC_folder);end % create folder if non-existant
	if paths.in_jack
	    AC_file = fullfile(paths.AC_folder,sprintf('Activation_%s_%s.mat',paths.results_jack_name,con_name));
	elseif isfield(paths,'results_sub_folder') && ~isempty(paths.results_sub_folder)
	    AC_file = fullfile(paths.AC_folder,sprintf('Activation_%s_%s_%s.mat',paths.results_name,paths.results_sub_folder,con_name));
	else
	    AC_file = fullfile(paths.AC_folder,sprintf('Activation_%s_Global_%s.mat',paths.results_name,con_name));
	end
	save(AC_file,'n_vox_mask','n_vox_maskP','meanT','meanTP','maxT','peak_threshold','Cluster_threshold','extent_threshold_vox','desired_p_val')
end

