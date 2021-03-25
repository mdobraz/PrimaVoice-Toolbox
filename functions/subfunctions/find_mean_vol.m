function ref_volume = find_mean_vol(Q,selected_vols,steadiest_period)

if ~exist('selected_vols','var') || isempty(selected_vols)
    selected_vols = true(size(Q,1),1);
end

if ~exist('steadiest_period','var') || isempty(steadiest_period)
    steadiest_period = true(size(Q,1),1);
end

mean_selected_trans = mean(Q(selected_vols,1:3)); % mean translations along the selected volumes
mean_selected_rot = mean(Q(selected_vols,4:6)); % mean rotations along the selected volumes
dist2mean_trans = sum(abs(Q(:,1:3) - repmat(mean_selected_trans,size(Q,1),1)),2); % distance of each volume to the mean translations
dist2mean_rot = sum(abs(Q(:,4:6) - repmat(mean_selected_rot,size(Q,1),1)),2); % distance of each volume to the mean rotations

dist2mean_trans(~selected_vols | ~steadiest_period) = inf; % exclude unselected volumes
dist2mean_rot(~selected_vols | ~steadiest_period) = inf;

[~,I_trans] = sort(dist2mean_trans);
[~,I_rot] = sort(dist2mean_rot);

rank_trans = zeros(length(dist2mean_trans),1);
rank_rot = zeros(length(dist2mean_trans),1);
for i = 1:length(dist2mean_trans)
    rank_trans(I_trans(i)) = i;
    rank_rot(I_rot(i)) = i;
end

[~,ref_volume] = min((rank_trans .* 2) + rank_rot); % most representative volume of the minimum movement period (put more weight on translations)