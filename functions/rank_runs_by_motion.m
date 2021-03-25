trans_rots = nan(n_total_runs,1);
nrejected = nan(n_total_runs,1);
session = nan(n_total_runs,1);
run = nan(n_total_runs,1);
labels = cell(n_total_runs,1);
run_count = 0;
for s = 1:numel(SR)
    runs = SR(s).runs;
    for r = 1:length(runs)
        run_count = run_count + 1;
        run(run_count) = runs(r);
        session(run_count) = SR(s).session;
        labels{run_count} = sprintf('%i-%i',session(run_count),run(run_count));
        % Filename
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);
        if strcmp(ext,'.gz'); [~,bold_name,ext] = fileparts(bold_name); end
        
        % Load RealignRef_file
        RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
        load(RealignRef_file)
        
        trans_rots(run_count) = norm([mcRMS.trans mcRMS.rot]); % norm of trans & rot
        nrejected(run_count) = sum(~selected_vols) / length(selected_vols);
        
    end
end


% trans_sorted = sort(trans,'ascend');
% [~,ranks_trans] = ismember(trans,trans_sorted); % ranks with equal ranks for equal results

% rot_sorted = sort(rots,'ascend');
% [~,ranks_rot] = ismember(rots,rot_sorted); % ranks with equal ranks for equal results

% classement = mean([ranks_trans ranks_rot],2);
% [~,MC_ranks] = sort(classement,'ascend');
% [~,classement] = sort(MC_ranks,'ascend');

TR_norm = nan(run_count,1);
for i = 1:run_count
    TR_norm(i) = norm([trans_rots(i) nrejected(i)]); % use norm instead of mean
end

[~,MC_ranks] = sort(TR_norm,'ascend');
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
xl = 'Translations, rotations norm';
yl = 'Ratio of rejected volumes';
titre = 'Runs motion RMS';

max_rank = round(length(MC_ranks) / 2);
plot_dyn_rank(ax,trans_rots,nrejected,MC_ranks,max_rank,xl,yl,titre)

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
                 'Callback','max_rank=round(rank_sl.Value); rank_sl_tx.String=max_rank;plot_dyn_rank(ax,trans_rots,nrejected,MC_ranks,max_rank,xl,yl,titre)');

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

fig_prefix = fullfile(paths.analysis,['MotionCorr_Run_selection_' num2str(max_rank)]);
figurewrite(fig_prefix,[],0,paths.analysis); % the 0 is to force eps figure
delete(rfig)
            

% construct matrix of selected runs
reduced_I = MC_ranks(1:max_rank);
reduced_sessions = session(reduced_I);
reduced_runs = run(reduced_I);
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


