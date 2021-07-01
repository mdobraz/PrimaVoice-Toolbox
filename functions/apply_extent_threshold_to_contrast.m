con = 'macaque_vs_all';
thres_file = fullfile(paths.results_multi,sprintf('sub-%s_res-%s_%s_threshold.mat',paths.subject,paths.results_name,con));
t = load(thres_file);

in_file = fullfile(paths.results_multi,sprintf('sub-%s_res-%s_%s_t.nii.gz',paths.subject,paths.results_name,con));
out_file = fullfile(paths.results_multi,sprintf('sub-%s_res-%s_%s_t_extent-thresholded.nii',paths.subject,paths.results_name,con));

if strcmp(paths.contrast_agent,'MION')
	isMION = 1;
else
	isMION = 0;
end

out_file = apply_extent_threshold(in_file,out_file,t.peak_threshold(1),t.Cluster_threshold,t.extent_threshold_vox,isMION);


% apply_extent_threshold_to_contrast
% out_file = xfm_func2anat(out_file,paths,coreg_params)
% system(['gzip ' out_file]);