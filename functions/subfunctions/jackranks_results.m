con_name = 'sound_vs_silence';
name_base = sprintf('Activation_%s',paths.results_name);

jack_files = dir(fullfile(paths.AC_folder,sprintf('%s_Jackknife_rank-*%s.mat',name_base,con_name)));

rank_extent = nan(numel(jack_files),1);
rank_Pextent = nan(numel(jack_files),1);
rank_Tmean = nan(numel(jack_files),1);
rank_TPmean = nan(numel(jack_files),1);
rank_Tmax = nan(numel(jack_files),1);
rank_ranks = nan(numel(jack_files),1);
rank_CT = nan(numel(jack_files),1);
rank_ETV = nan(numel(jack_files),1);


for f = 1:numel(jack_files)
    fname = jack_files(f).name;
    load(fullfile(paths.AC_folder,fname))
    rank_ranks(f) = str2double(fname(strfind(fname,'_rank-')+6:strfind(fname,'_rank-')+8));
    rank_extent(f) = n_vox_mask;
    rank_Pextent(f) = n_vox_maskP;
    rank_Tmean(f) = meanT - Cluster_threshold;
    rank_TPmean(f) = meanTP - peak_threshold;
    rank_Tmax(f) = maxT - peak_threshold;
    rank_CT(f) = Cluster_threshold;
    rank_ETV(f) = extent_threshold_vox;
end


return



linsp = (1:length(rank_Pextent))';
Pextent_sorted = sort(rank_Pextent,'ascend');
[~,ranks_Pextent] = ismember(rank_Pextent,Pextent_sorted); % ranks with equal ranks for equal results

Tmax_sorted = sort(rank_Tmax,'ascend');
[~,ranks_Tmax] = ismember(rank_Tmax,Tmax_sorted); % ranks with equal ranks for equal results

classement = mean([ranks_Pextent ranks_Tmax],2);
[~,jack_ranks] = sort(classement,'ascend');
[~,classement] = sort(jack_ranks,'ascend');



% rank_labels(I)

