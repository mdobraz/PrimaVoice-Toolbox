function [ref_volume, selected_vols, Q, i1, i2] = find_refvol_SPM(Scans,flags,mvt_params,paths,mag_file,in_file)

%% Get realignment estimates relative to the first image
flags.rtm = 0; % no registration to mean
P = spm_realign(Scans,flags);

%% Transform to a readable matrix
n = length(P);
Q = zeros(n,6);
for j=1:n
    qq     = spm_imatrix(P(j).mat/P(1).mat);
    Q(j,:) = qq(1:6);
end
Q(:,4:6) = rad2deg(Q(:,4:6)); % radians to degrees


%% volumes with minimal translation
trans_Q_diff = abs(diff(Q(:,1:3))); % absolute difference of the mvt estimate
min_translation_vols = trans_Q_diff(:,1) < mvt_params.max_trans_x & trans_Q_diff(:,2) < mvt_params.max_trans_y & trans_Q_diff(:,3) < mvt_params.max_trans_z; % indices of volumes with minimal translation

%% volumes with minimal rotation
rot_Q_diff = abs(diff(Q(:,4:6))); % absolute difference of the rotation estimate
min_rotation_vols = rot_Q_diff(:,1) < mvt_params.max_rot_pitch & rot_Q_diff(:,2) < mvt_params.max_rot_roll & rot_Q_diff(:,3) < mvt_params.max_rot_yaw; % indices of volumes with minimal rotation

%% selected volumes is the intersection between min rot and trans volumes
selected_vols = min_rotation_vols & min_translation_vols;
selected_vols = [selected_vols(1);selected_vols]; % add the first volume

%% remove groups of consecutive ones or zeros in X, whose size is below 'min_n_consec_ones' or 'min_n_consec_zeros'
min_n_consec_ones = ceil(mvt_params.min_dur_steady/(mvt_params.TR)); % in volumes, number of consecutive ones (steady)
min_n_consec_zeros = ceil(mvt_params.min_dur_mvt/(mvt_params.TR)); % in volumes, number of consecutive zeros (movement)
selected_vols = clean_binary(selected_vols,min_n_consec_ones,min_n_consec_zeros);

%% Refine selection by removing volumes to far from the longuest & steadiest period
if any(selected_vols)
    group_sizes = get_group_sizes(selected_vols);
    for g = 1:size(group_sizes,1)
        if group_sizes(g,2)
            if g > 1
                i1 = sum(group_sizes(1:g-1,1)) + 1;
                i2 = i1 + group_sizes(g,1) - 1;
            else
                i1 = 2;
                i2 = i1 + group_sizes(g,1) - 2;
            end
            i1 = i1 - 1; % because working with Q_diffs
            i2 = i2 - 1;
            group_sizes(g,3) = sum(mean(trans_Q_diff(i1:i2,:),1));
            group_sizes(g,4) = sum(mean(rot_Q_diff(i1:i2,:),1));
        end
    end

    [~,best_group] = max(prod(group_sizes(:,1:2),2) ./ prod(group_sizes(:,3:4),2)); % group of the longuest & steadiest period
    if best_group > 1
        i1 = sum(group_sizes(1:best_group-1,1)) + 1;
    else
        i1 = 1;
    end
    i2 = i1 + group_sizes(best_group,1) - 1;

    steadiest_period = false(length(selected_vols),1);
    steadiest_period(i1:i2) = true;

    mean_selected_trans = mean(Q(steadiest_period,1:3)); % mean translations along the selected volumes
    mean_selected_rot = mean(Q(steadiest_period,4:6)); % mean rotations along the selected volumes
    dist2mean_trans = max(abs(Q(:,1:3) - repmat(mean_selected_trans,length(Q),1)),[],2); % distance of each volume to the mean translations
    dist2mean_rot = max(abs(Q(:,4:6) - repmat(mean_selected_rot,length(Q),1)),[],2); % distance of each volume to the mean rotations

    selected_vols = selected_vols & (dist2mean_trans < (max([mvt_params.max_trans_x mvt_params.max_trans_y mvt_params.max_trans_z]) * mvt_params.trans_fact));
    selected_vols = selected_vols & (dist2mean_rot < (max([mvt_params.max_rot_pitch mvt_params.max_rot_roll mvt_params.max_rot_yaw]) * mvt_params.rot_fact));

    %% clean again after selection by absolute coords
    selected_vols = clean_binary(selected_vols,min_n_consec_ones,min_n_consec_zeros);

    %% Erode selected_vols (remove volumes at the frontier of the movement)
    for i = 1:mvt_params.after_mvt_lag % erode by a factor of 'after_mvt_lag'
        selected_vols = erode_binary(selected_vols,'before'); % before steady period
    end

    %% clean again after erosion
    selected_vols = clean_binary(selected_vols,min_n_consec_ones,min_n_consec_zeros);

    %% find representative volume
    if exist('mag_file','var') && mvt_params.sel_from_fmap
        fprintf('Selecting volume the closet to the fieldmap...')
        normmis = zeros(numel(selected_vols),1);
        roi_out_file = fullfile(paths.tmp_nii,'one_vol');
        for i = 1:numel(selected_vols)
            if steadiest_period(i)
                system(sprintf('%sfslroi %s %s %i 1',paths.FSL_prefix,in_file,roi_out_file,i-1)); % fslroi indexing starts at 0, hence the '-1' here
                [~,cost] = system(sprintf('%sflirt -in %s -ref %s -schedule ./functions/measurecost1.sch -cost normmi | head -1 | cut -f1 -d'' ''',paths.FSL_prefix,roi_out_file,mag_file));
                normmis(i) = str2double(cost);
            end
        end
        [~,ref_volume] = min(normmis);
        fprintf(' done.\n')
    else
        ref_volume = find_mean_vol(Q,selected_vols,steadiest_period); % most representative volume of the minimum movement period
    end
else
    ref_volume = nan;
    fprintf('No volume was selected in this run !!!')
end





    

















