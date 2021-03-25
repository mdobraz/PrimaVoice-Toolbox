con_name = 'sound_vs_silence';
name_base = sprintf('Activation_%s',paths.results_name);

jack_files = dir(fullfile(paths.AC_folder,sprintf('%s_Jackknife_ses-*%s.mat',name_base,con_name)));
jack_files(end+1) = jack_files(1);
jack_files(end).name = sprintf('%s_Global_%s.mat',name_base,con_name);

numjack = numel(jack_files);

extent = nan(numjack,1);
Pextent = nan(numjack,1);
Tmean = nan(numjack,1);
TPmean = nan(numjack,1);
Tmax = nan(numjack,1);
sessions = nan(numjack,1);
runs = nan(numjack,1);
labels = cell(numjack,1);

% sel_vols_ratio = nan(numel(files),1);
% ok_blocks_ratio= nan(numel(files),1);



for f = 1:numjack
    fname = jack_files(f).name;
    load(fullfile(paths.AC_folder,fname))
    if strcmp(fname,[name_base '_Global_' con_name '.mat'])
        all_runs = f;
        sessions(f) = inf;
        runs(f) = inf;
        labels{f} = 'All runs';
%         sel_vols_ratio(f) = inf;
%         ok_blocks_ratio(f) = inf;
    else
        sessions(f) = str2double(fname(strfind(fname,'_ses-')+5:strfind(fname,'_ses-')+6));
        runs(f) = str2double(fname(strfind(fname,'-run-')+5:strfind(fname,'-run-')+6));
        labels{f} = sprintf('%i-%i',sessions(f),runs(f));
%         check_file = sprintf('%s/check_%s_session%i_bold%i.mat',paths.cond_mat,paths.subject,sessions(f),runs(f));
%         load(check_file)
%         ok_blocks_ratio(f) = sum(n_ok) / sum(n_pres);
%         sel_vols_ratio(f) = sum(selected_vols) / length(selected_vols);
    end
    extent(f) = n_vox_mask;
    Pextent(f) = n_vox_maskP;
    Tmean(f) = meanT - Cluster_threshold;
    TPmean(f) = meanTP - peak_threshold;
    Tmax(f) = maxT - peak_threshold;
end

% sel_vols_ratio(sel_vols_ratio == inf) = mean(sel_vols_ratio(sel_vols_ratio ~= inf));
% ok_blocks_ratio(ok_blocks_ratio == inf) = mean(ok_blocks_ratio(ok_blocks_ratio ~= inf));

% figure
% plot(sel_vols_ratio(Pextent > Pextent(runs == inf)),ok_blocks_ratio(Pextent > Pextent(runs == inf)),'ro','MarkerFaceColor','r')
% hold on
% grid on
% plot(sel_vols_ratio(Pextent <= Pextent(runs == inf)),ok_blocks_ratio(Pextent <= Pextent(runs == inf)),'ko','MarkerFaceColor','k')
% xlabel('selected vols ratio')
% ylabel('ok blocks ratio')
% hold off

% figure
% x = 1:length(Pextent);
% plot(x,Pextent,'o','MarkerFaceColor','k')
% hold on
% grid on
% set(gca,'xTick',x)
% set(gca,'xTickLabel',labels)
% xrange = get(gca,'xlim');
% plot(xrange,[Pextent(runs==inf) Pextent(runs==inf)])
% ylabel('Peak Extent')


% figure
% x = 1:length(ok_blocks_ratio);
% plot(x,ok_blocks_ratio,'o','MarkerFaceColor','k')
% hold on
% grid on
% set(gca,'xTick',x)
% set(gca,'xTickLabel',labels)
% xrange = get(gca,'xlim');
% plot(xrange,[ok_blocks_ratio(runs==inf) ok_blocks_ratio(runs==inf)])
% ylabel('Ok blocks ratio')


% linsp = (1:length(Pextent))';
% Pextent_sorted = sort(Pextent,'ascend');
% [~,ranks_Pextent] = ismember(Pextent,Pextent_sorted); % ranks with equal ranks for equal results

% Tmax_sorted = sort(Tmax,'ascend');
% [~,ranks_Tmax] = ismember(Tmax,Tmax_sorted); % ranks with equal ranks for equal results

% % ok_blocks_ratio_sorted = sort(ok_blocks_ratio,'descend');
% % [~,ranks_ok_ratio] = ismember(ok_blocks_ratio,ok_blocks_ratio_sorted); % ranks with equal ranks for equal results

% % classement = mean([ranks_Pextent ranks_ok_ratio ranks_Tmax],2);
% classement = mean([ranks_Pextent ranks_Tmax],2);
% [~,jack_ranks] = sort(classement,'ascend');
% [~,classement] = sort(jack_ranks,'ascend');
% classement(runs == Inf) = []; % remove all runs in ranking
% % labels(I)





TR_norm = nan(numjack,1);
for i = 1:numjack
    TR_norm(i) = norm([Pextent(i) Tmax(i)]); % use norm instead of mean
end

[~,jack_ranks] = sort(TR_norm,'ascend');
[~,classement] = sort(jack_ranks,'ascend');
classement(runs == Inf) = []; % remove all runs in ranking






run_count = 0;
all_ranks = cell(numel(SR),1);
for s = 1:numel(SR)
    SRruns = SR(s).runs;
    all_ranks{s} = nan(1,length(SRruns));
    for r = 1:length(SRruns)
        run_count = run_count + 1;
        all_ranks{s}(r) = classement(run_count);
    end
end

linsp = (1:length(Pextent))';