function [jack_sessions,jack_runs,rem_session,rem_run] = get_sessions_runs_jackknife(SR,run_index)

jack_sessions = [SR.session];
nsessions = numel(SR);
jack_runs = cell(nsessions,1);

for s = 1:nsessions
    jack_runs{s} = SR(s).runs;
end

ri = 0;
stop = 0;
for s = 1:nsessions
    for r = 1:length(jack_runs{s})
        ri = ri + 1;
        if ri == run_index
            rem_session = jack_sessions(s);
            these_runs = jack_runs{s};
            rem_run = these_runs(r);
            these_runs(r) = [];
            if isempty(these_runs)
                jack_sessions(s) = [];
                jack_runs(s) = [];
            else
                jack_runs{s} = these_runs;
            end
            stop = 1;
        end
        if stop;break;end
    end
    if stop;break;end
end


fprintf('\n\n### Session(s):\n')
fprintf('%02.0f ',jack_sessions)

fprintf('\n### Run(s):\n')
for s = 1:length(jack_runs)
    fprintf('%02.0f ',jack_runs{s})
    fprintf('\n')
end

fprintf('\nRemoved run: ses-%02.0f run-%02.0f\n\n',rem_session,rem_run)


