fprintf('\nInter session cost matrix...\n\n')
% Retrieve files
avg_vols = cell(n_total_runs,1);
cost_sessions = nan(n_total_runs,1);
cost_runs = nan(n_total_runs,1);
nrun = 0;
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        nrun = nrun + 1;
        cost_sessions(nrun) = session;
        cost_runs(nrun) = run;
        bold_file = SR(s).filename{r};
        [~,bold_name,ext] = fileparts(bold_file);
        if strcmp(ext,'.gz'); [~,bold_name,ext] = fileparts(bold_name); end
        
        
        
        name_base = [reslice_flags.prefix 'preMC_' bold_name];
        if realign_FSL.use_avg
            avg_vols{nrun} = fullfile(paths.resliced,sprintf('flirt_%s_avg',name_base));
        else
            RealignRef_file = fullfile(paths.realign,[bold_name '_RealignRef.mat']);
            load(RealignRef_file)
            avg_vols{nrun} = fullfile(paths.resliced,sprintf('flirt_%s_%i',name_base,ref_volume));
        end
    end
end

% Create cost matrix
cost_matrix = nan(n_total_runs,n_total_runs);
cost_avg = nan(1,n_total_runs);
for i = 1:n_total_runs % i is the ref
    fprintf('run %i / %i\n',i,n_total_runs)
    for j = 1:n_total_runs
        [~,cost] = system(sprintf('%sflirt -in %s -ref %s -schedule ./functions/measurecost1.sch -cost %s | head -1 | cut -f1 -d'' ''',paths.FSL_prefix,avg_vols{j},avg_vols{i},realign_FSL.flirt_cost));
        cost_matrix(i,j) = str2double(cost);
    end
    [~,cost] = system(sprintf('%sflirt -in %s -ref %s -schedule ./functions/measurecost1.sch -cost %s | head -1 | cut -f1 -d'' ''',paths.FSL_prefix,avg_vols{i},paths.average,realign_FSL.flirt_cost));
        cost_avg(i) = str2double(cost);
end

cost_matrix(logical(eye(n_total_runs))) = nan;


% Save figure
figureprep([0 0 800 800]);
% figure('position',[100 100 600 800])
sp1 = subplot(2,1,1);
imagesc(cost_matrix)
sp2 = subplot(2,1,2);
imagesc(cost_avg)
sp1p = sp1.Position;
sp2p = sp2.Position;
sp1.Position = [sp1p(1) (sp2p(2) + 0.15) sp1p(3) 0.7];
sp2.Position = [sp2p(1) sp2p(2) sp1p(3) 0.05];
fig_prefix = fullfile(paths.analysis,'inter_session_cost_matrix');
figurewrite(fig_prefix,[],0,paths.analysis); % the 0 is to force eps figure

% Compute ranking
mean1 = nanmean(cost_matrix,1);
mean2 = nanmean(cost_matrix,2)';
gmean = (mean1 + mean2 + (2 * cost_avg)) / 4;

[~,costrank] = sort(gmean);

% save
save(paths.costmatrix,'cost_matrix','costrank','cost_sessions','cost_runs','cost_avg')

% Write merged averages
command = sprintf('%sfslmerge -t %s',paths.FSL_prefix,fullfile(paths.resliced,'Avg_Verif'));
for i = 1:n_total_runs
    command = [command ' ' avg_vols{i}];
end
system(command);








