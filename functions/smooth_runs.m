%% Smooth all sessions/runs
%%%%%% Run 'set_parameters' & 'realign_runs' beforehand

%% SCRIPT
if ~exist(paths.smoothed,'dir');mkdir(paths.smoothed);end % create folder if non-existant
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        fprintf('\nSmoothing session %i, run %i:\n',session,run)
        [~,name,ext] = fileparts(SR(s).filename{r});
        resliced_file = fullfile(paths.resliced,[reslice_flags.prefix name ext]);
        smoothed_file = fullfile(paths.smoothed,[smooth_params.prefix reslice_flags.prefix name ext]);
        if mean(mvt_params.voxel_size == smooth_params.fwhm) == 1 % no smoothing necessary / copy file
            copyfile(resliced_file,smoothed_file);
        else % smooth runs
            for vol = 1:scan_info(s,r).NumberOfVolumesInFile
                fprintf('.')
                orig = sprintf('%s,%i',resliced_file,vol);
                smoothed = sprintf('%s,%i',smoothed_file,vol);
                spm_smooth(orig,smoothed,smooth_params.fwhm) % perform smoothing
                if ~mod(vol,100)
                    fprintf('\n')
                end
            end
        end
    end
end



