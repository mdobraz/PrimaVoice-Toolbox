SR_ME = struct_sess_run(BIDS,'all','all',paths.subject,[],'epi','anat');


n_ME = 0;
for s = 1:numel(SR_ME)
    session = SR_ME(s).session;
    runs = SR_ME(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        n_ME = n_ME + 1;
        verif_files{n_ME} = compute_T2s(SR_ME(s).filename{r},session,run,paths,BIDS);
    end
end


merge_in_files = sprintf('%s ',verif_files{:});
merge_out_file = fullfile(paths.ME,'merged_brains.nii.gz');
system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,merge_out_file,merge_in_files)); % merge anat runs

% ME_data = spm_BIDS(BIDS,'data','sub',paths.subject,'type','epi','modality','anat');