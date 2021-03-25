
mask_contrast = RDM_params.contrast;
peak_threshold = RDM_params.peak_threshold;
mask_with_grey = RDM_params.mask_with_grey;
min_clust_size = RDM_params.min_clust_size;
n_voxels = RDM_params.n_voxels;
n_peaks = RDM_params.n_peaks;
ScanLog = ScanLog_ref;
pdistance = RDM_params.distance;
colmap = RDM_params.colormap;


if isfield(paths,'subjects')
    fixed_random = RDM_params.fixed_random;
    FRstr = [fixed_random '_'];
else
    FRstr = '';
end

if n_voxels > min_clust_size
    min_clust_size = n_voxels;
    warning('Needed neighbors is lower than minimal cluster extent. min_clust_size is now %d.\n', n_voxels);
end


%% load RDM file
thr_str = sprintf('%g',peak_threshold);
thr_str = strrep(thr_str,'.','p');
peaks_str = sprintf('%g',n_peaks);
peaks_str = strrep(peaks_str,'.','p');

if mask_with_grey
    grey_str = '_grey-masked';
else
    grey_str = '';
end
if isfield(paths,'CV_results_sub_folder') && ~isempty(paths.CV_results_sub_folder)
    CV_str = ['_CV-' paths.CV_results_sub_folder];
    isCV = 1;
else
    CV_str = '';
    isCV = 0;
end
if ~isfield(paths,'L') || isempty(paths.L)
    L_str = sprintf('_L%i',numel(fieldnames(ScanLog.fMRISTAT)));
else
    L_str = sprintf('_L%i',paths.L);
end

RDM_name = sprintf('%s%s_brain_RDM_%s_thres-%s_clust-size-%i_n-vox-%i_n-peaks-%s%s%s%s',FRstr,pdistance,mask_contrast,thr_str,min_clust_size,n_voxels,peaks_str,grey_str,CV_str,L_str);

RDMs_dir = fullfile(paths.results_multi,'brain_RDMs');
if ~exist(RDMs_dir,'dir');mkdir(RDMs_dir);end % create tmp folder if non-existant
RDM_dir = fullfile(RDMs_dir,RDM_name);
if ~exist(RDM_dir,'dir');mkdir(RDM_dir);end % create tmp folder if non-existant

if isfield(paths,'subjects')
    prefix = ['second_level_' paths.analysis_name '_' RDM_name];
else
    prefix = ['sub-' paths.subject '_res-' paths.results_name '_' RDM_name];
end

RDM_file = fullfile(RDM_dir,[prefix '.mat']);

if isfield(paths,'subject')
    subject = ['sub-' paths.subject];
else
    subject = ['second_level_' paths.analysis_name];
end




