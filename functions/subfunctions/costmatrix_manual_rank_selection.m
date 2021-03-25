CR = load(paths.costmatrix);
sessions = CR.cost_sessions;
runs = CR.cost_runs;
jack_ranks = CR.costrank;

% load activations
con_name = 'sound_vs_silence';
name_base = sprintf('Activation_%s',paths.results_name);

cost_files = dir(fullfile(paths.AC_folder,sprintf('%s_Costmatrix_rank-*%s.mat',name_base,con_name)));

rank_extent = nan(numel(cost_files),1);
rank_Pextent = nan(numel(cost_files),1);
rank_Tmean = nan(numel(cost_files),1);
rank_TPmean = nan(numel(cost_files),1);
rank_Tmax = nan(numel(cost_files),1);
rank_ranks = nan(numel(cost_files),1);
rank_CT = nan(numel(cost_files),1);
rank_ETV = nan(numel(cost_files),1);
labels = cell(numel(cost_files),1);


for f = 1:numel(cost_files)
    fname = cost_files(f).name;
    load(fullfile(paths.AC_folder,fname))
    rank_ranks(f) = str2double(fname(strfind(fname,'_rank-')+6:strfind(fname,'_rank-')+8));
    rank_extent(f) = n_vox_mask;
    rank_Pextent(f) = n_vox_maskP;
    rank_Tmean(f) = meanT - Cluster_threshold;
    rank_TPmean(f) = meanTP - peak_threshold;
    rank_Tmax(f) = maxT - peak_threshold;
    rank_CT(f) = Cluster_threshold;
    rank_ETV(f) = extent_threshold_vox;
    labels{f} = sprintf('%i-%i',sessions(f),runs(f));
end

linsp = (1:length(rank_Pextent))';


g = groot;
if size(g.MonitorPositions,1) > 1
    sx = g.MonitorPositions(2,1) / abs(g.MonitorPositions(2,1));
else
    sx = 0;
end
figure('Units','Normalized','Position',[sx 0 1 1])
flinsp = flip(linsp);

hold on
for i = 1:length(linsp)
    text(0.2,flinsp(i),labels(jack_ranks(i)))
end
axis([0 2 0 length(linsp)+1])
set(gca,'yTick',linsp)
set(gca,'yTickLabel',flinsp)


figure('Units','Normalized','Position',[0.05 0.3 0.4 0.4])
plot(rank_Pextent,'-o')
xlabel('Ranks included')
ylabel('Sound vs Silence spatial extent')
grid on

figure('Units','Normalized','Position',[0.5 0.3 0.4 0.4])
plot(rank_Tmax,'-o')
xlabel('Ranks included')
ylabel('Sound vs Silence max t-value')
grid on


max_rank = str2double(input('At which rank do you want to stop selection (included)? ','s'));

% construct matrix of selected runs
reduced_I = jack_ranks(1:max_rank);
reduced_sessions = sessions(reduced_I);
reduced_sessions = reduced_sessions(reduced_sessions ~= inf);
reduced_runs = runs(reduced_I);
reduced_runs = reduced_runs(reduced_runs ~= inf);
selected_sessions = unique(reduced_sessions);
selected_runs = cell(length(selected_sessions),1);

for i = 1:length(selected_sessions)
    sess = selected_sessions(i);
    selected_runs{i} = sort(reduced_runs(reduced_sessions == sess));
end

fprintf('\nAll selected runs:\n')
for i = 1:length(selected_runs)
    fprintf('%i\t',selected_sessions(i))
    fprintf('%i ',selected_runs{i})
    fprintf('\n');
end


% Save
filename = sprintf('%s/Costmatrix_Run_selection.mat',paths.analysis);
save(filename,'selected_sessions','selected_runs')
fprintf('\n%s ....... saved\n\n',filename)





% create several lists of runs with equal ranks
n_groups = str2double(input('In how many groups do you want to separate the runs? ','s'));
selection = select_equal_ranks(max_rank,n_groups,10000);

for g = 1:n_groups
    % construct matrix of selected runs
    reduced_Ig = reduced_I(selection(:,g));
    reduced_sessions = sessions(reduced_Ig);
    reduced_sessions = reduced_sessions(reduced_sessions ~= inf);
    reduced_runs = runs(reduced_Ig);
    reduced_runs = reduced_runs(reduced_runs ~= inf);
    selected_sessions = unique(reduced_sessions);
    selected_runs = cell(length(selected_sessions),1);

    for i = 1:length(selected_sessions)
        sess = selected_sessions(i);
        selected_runs{i} = sort(reduced_runs(reduced_sessions == sess));
    end

    fprintf('\nSelected runs group %i:\n',g)
    for i = 1:length(selected_runs)
        fprintf('%i\t',selected_sessions(i))
        fprintf('%i ',selected_runs{i})
        fprintf('\n');
    end


    % Save
    filename = sprintf('%s/Costmatrix_Run_selection_%i_folds_group%i.mat',paths.analysis,n_groups,g);
    save(filename,'selected_sessions','selected_runs')
    fprintf('\n%s ....... saved\n\n',filename)
end
