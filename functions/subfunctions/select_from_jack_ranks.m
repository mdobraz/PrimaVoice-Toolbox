function [selected_sessions,selected_runs,reduced_I,reduced_sessions,reduced_runs] = select_from_jack_ranks(sessions,runs,jack_ranks,max_rank)

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