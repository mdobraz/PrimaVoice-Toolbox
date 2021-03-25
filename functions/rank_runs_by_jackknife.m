jackknife_results
jackranks_results

g = groot;
if size(g.MonitorPositions,1) > 1
    sx = g.MonitorPositions(2,1) / abs(g.MonitorPositions(2,1));
else
    sx = 0;
end
gf1 = figure('Units','Normalized','Position',[sx 0 1 1]);
flinsp = flip(linsp);

hold on
for i = 1:length(linsp)
%     text(1,flinsp(i),num2str(linsp(i)))
%     if (Pextent(I(i)) > Pextent(runs == inf)) && (Tmax(I(i)) > Tmax(runs == inf)) && (ok_blocks_ratio(I(i)) < ok_blocks_ratio(runs == inf))
    if (Pextent(jack_ranks(i)) > Pextent(runs == inf)) && (Tmax(jack_ranks(i)) > Tmax(runs == inf))
        color = 'r';
    elseif (Pextent(jack_ranks(i)) < Pextent(runs == inf)) && (Tmax(jack_ranks(i)) < Tmax(runs == inf))
        color = 'g';
    else
        color = 'k';
    end
    text(0.2,flinsp(i),labels(jack_ranks(i)),'Color',color)
end
axis([0 2 0 length(linsp)+1])
set(gca,'yTick',linsp)
set(gca,'yTickLabel',flinsp)


gf2 = figure('Units','Normalized','Position',[0.05 0.3 0.4 0.4]);
plot(rank_Pextent,'-o')
xlabel('Ranks included')
ylabel('Sound vs Silence spatial extent')
grid on

gf3 = figure('Units','Normalized','Position',[0.5 0.3 0.4 0.4]);
plot(rank_Tmax,'-o')
xlabel('Ranks included')
ylabel('Sound vs Silence max t-value')
grid on


max_rank = str2double(input('At which rank do you want to stop selection (included)? ','s'));

close(gf1);close(gf2);close(gf3)

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


figure
plot(Pextent(jack_ranks(1:max_rank)),Tmax(jack_ranks(1:max_rank)),'ko','MarkerFaceColor','b')
hold on
grid on
plot(Pextent(jack_ranks(max_rank+2:end)),Tmax(jack_ranks(max_rank+2:end)),'ko','MarkerFaceColor','r')
plot(Pextent(end),Tmax(end),'ko','MarkerFaceColor','g','MarkerSize',10)
xlabel('Spatial extent (voxels)')
ylabel('Max t-value')
title('Press enter to continue (do not close this figure)')
hold off

fig_prefix = fullfile(paths.analysis,['Jackknife_Run_selection_' num2str(max_rank)]);
input('### Press enter to close figure and continue script')
figurewrite(fig_prefix,[],0,paths.analysis); % the 0 is to force eps figure


% Save
filename = [fig_prefix '.mat'];
save(filename,'selected_sessions','selected_runs','all_ranks','classement','max_rank')
fprintf('\n%s ....... saved\n\n',filename)





% create several lists of runs with equal ranks
n_groups = str2double(input('In how many groups do you want to separate the runs? ','s'));
values = Pextent(reduced_I) .* Tmax(reduced_I);
selection = select_equal_ranks(values,n_groups,10000);

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
    filename = [fig_prefix '.mat'];
    filename = sprintf('%s_%i_folds_group%i.mat',fig_prefix,n_groups,g);
    save(filename,'selected_sessions','selected_runs')
    fprintf('\n%s ....... saved\n\n',filename)
end

