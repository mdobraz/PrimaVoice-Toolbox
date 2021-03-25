before = nan(n_total_runs,1);
after = nan(n_total_runs,1);
sessionsR = nan(n_total_runs,1);
runsR = nan(n_total_runs,1);
labels = cell(n_total_runs,1);
run_count = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    % get T2s
    times = nan(length(runs),1);
    mean_R = nan(length(runs),1);
    for r = 1:length(runs)
        run = runs(r);
        T2s_file = fullfile(paths.ME,sprintf('sub-%s_ses-%02.0f_run-%02.0f_T2s.mat',paths.subject,session,run));
		T2s = load(T2s_file);
        mean_R(r) = T2s.mean_R;
		times(r) = datenum(T2s.scan_time);
    end
    
    if length(times) > 2
        b = robustfit(times,mean_R);
    else
        b = regress(mean_R,[ones(length(times),1) times]);
    end
    
    for r = 1:length(runs)
        run_count = run_count + 1;
        runsR(run_count) = runs(r);
        sessionsR(run_count) = SR(s).session;
        labels{run_count} = sprintf('%i-%i',sessionsR(run_count),runsR(run_count));
        % Bold time
        bold_AT = scan_info(s,r).AcquisitionDateTime;
        bold_time = datenum(datetime(bold_AT,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS'));
        
        ibef = find(times<bold_time,1,'last');
        iaft = find(times>bold_time,1,'first');
        
        if isempty(ibef) % deal with forgotten scans...
            ibef = bold_time - 0.001;
            after(run_count) = b(1)+b(2)*ibef;
        else
            before(run_count) = mean_R(ibef);
        end
        if isempty(iaft) % deal with forgotten scans...
            durrun = scan_log(s,r).RM1_trigs(end);
            iaft = addtodate(bold_time,ceil(durrun), 'second');
            after(run_count) = b(1)+b(2)*iaft;
        else
            after(run_count) = mean_R(iaft);
        end
    end
end


% before_sorted = sort(before,'descend');
% [~,ranks_before] = ismember(before,before_sorted); % ranks with equal ranks for equal results

% after_sorted = sort(after,'descend');
% [~,ranks_after] = ismember(after,after_sorted); % ranks with equal ranks for equal results

% classement = mean([ranks_before ranks_after],2);
% [~,MC_ranks] = sort(classement,'ascend');
% [~,classement] = sort(MC_ranks,'ascend');


BA_norm = nan(run_count,1);
for i = 1:run_count
    BA_norm(i) = norm([before(i) after(i)]); % use norm instead of mean
end

[~,MC_ranks] = sort(BA_norm,'descend');
[~,classement] = sort(MC_ranks,'ascend');



run_count = 0;
all_ranks = cell(numel(SR),1);
for s = 1:numel(SR)
    runs = SR(s).runs;
    all_ranks{s} = nan(1,length(runs));
    for r = 1:length(runs)
        run_count = run_count + 1;
        all_ranks{s}(r) = classement(run_count);
    end
end

g = groot;
if size(g.MonitorPositions,1) > 1
    sx = g.MonitorPositions(2,1) / abs(g.MonitorPositions(2,1));
else
    sx = 0;
end

rfig = figure('Units','Normalized','Position',[sx 0 1 1]);
ax = axes('Parent',rfig);

% p = plot(trans,rots,'.','Parent',rax);
% text(trans(MC_ranks) + max(trans)/200,rots(MC_ranks) + max(rots)/200,num2cell(classement(MC_ranks)),'Parent',rax)
xl = 'R_2* before bold run';
yl = 'R_2* after bold run';
titre = 'Relaxation rates';

max_rank = round(length(MC_ranks) / 2);
plot_dyn_rank(ax,before,after,MC_ranks,max_rank,xl,yl,titre)

gui = figure('Toolbar','none',...
    'units','Normalized',...
    'Name','Selected the max rank',...
    'MenuBar','none',...
    'NumberTitle','off',...
    'Position',[0.4 0.4 0.2 0.1],...
    'resize','off');


rank_sl = uicontrol('style','slider',...
                 'parent',gui,...
                 'unit','normalized',...
                 'Value',max_rank,...
                 'Min',1,...
                 'Max',length(MC_ranks),...
                 'SliderStep',[1 / (length(MC_ranks) - 1) 5 / (length(MC_ranks) - 1)],...
                 'position',[0.1 0.5 0.8 0.2],...
                 'Callback','max_rank=round(rank_sl.Value); rank_sl_tx.String=max_rank;plot_dyn_rank(ax,before,after,MC_ranks,max_rank,xl,yl,titre)');

rank_sl_tx = uicontrol('style','text',...
                 'parent',gui,...
                 'units','normalized',...
                 'FontSize',(g.MonitorPositions(1,3) / 256) * 2,...
                 'HorizontalAlignment','center',...
                 'String',num2str(max_rank),...
                 'position',[0.4 0.2 0.2 0.2]);

rank_cl_tx = uicontrol('style','text',...
                 'parent',gui,...
                 'units','normalized',...
                 'FontSize',(g.MonitorPositions(1,3) / 256) * 1.5,...
                 'HorizontalAlignment','center',...
                 'String','Close this window when finished',...
                 'position',[0.1 0.75 0.8 0.2]);
uiwait(gui)

fig_prefix = fullfile(paths.subject_analysis,['R2s_Run_selection_' num2str(max_rank)]);
figurewrite(fig_prefix,[],0,paths.subject_analysis); % the 0 is to force eps figure
delete(rfig)
            

% construct matrix of selected runs
reduced_I = MC_ranks(1:max_rank);
reduced_sessions = sessionsR(reduced_I);
reduced_runs = runsR(reduced_I);
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
filename = [fig_prefix '.mat'];
save(filename,'selected_sessions','selected_runs','all_ranks','classement','max_rank')
fprintf('\nSelection saved in:\n%s\n\n',filename)


