function [RDMs,kRDMs,magnitudes] = plot_mean_brain_RDMs(RDM_files,opt)
% RDM_files: cell array of path to RDMs
% clust_nums: n_RDM_files x n_clust matrix of cluster numbers


if exist('opt','var')
    if isfield(opt,'clust_nums')
        clust_nums = opt.clust_nums;
    end
    if isfield(opt,'clust_names')
        clust_names = opt.clust_names;
    end
    if isfield(opt,'stim_names')
        stim_names = opt.stim_names;
    end
    if isfield(opt,'distances')
        distances = opt.distances;
    end
    if isfield(opt,'colors')
        def_colors = opt.colors;
    end
    if isfield(opt,'colormap')
        colmap = opt.colormap;
    else
        colmap = 'parula';
    end
    if isfield(opt,'dispFigs')
        dispFigs = opt.dispFigs;
    else
        dispFigs = 1;
    end
end






%% load RDM files and average them
for i = 1:numel(RDM_files)
    RDMs = load(RDM_files{i});

    if i == 1
        n_stims = size(RDMs.RDM,1);
        if exist('clust_nums','var')
            nclust = size(clust_nums,2);
        else
            nclust = length(RDMs.cluster);
            clust_nums = repmat(1:nclust,numel(RDM_files),1);
        end
        if ~exist('stim_names','var')
            stim_names = 1:n_stims;
        end
        RDM_avg = RDMs.RDM(:,:,clust_nums(i,:));

        medians = nan(numel(RDM_files),n_stims,nclust);
        means = nan(numel(RDM_files),n_stims,nclust);

        if ~exist('clust_names','var')
            clust_names = cell(nclust,1);
            for clust = 1:nclust
                clust_names{clust} = sprintf('Cluster %i',clust);
            end
        end
    else
        RDM_avg = RDM_avg + RDMs.RDM(:,:,clust_nums(i,:));
    end
    
    for c = 1:nclust
        kRDMs(c,i).RDM = RDMs.RDM(:,:,clust_nums(i,c));
        kRDMs(c,i).name = clust_names{c};
        kRDMs(c,i).color = [0 0 0];
    end

    for clust = 1:nclust
        medians(i,:,clust) = [RDMs.cluster(clust_nums(i,clust)).stim(:).medianT];
        means(i,:,clust) = [RDMs.cluster(clust_nums(i,clust)).stim(:).meanT];
    end

end
RDM_avg = RDM_avg ./ numel(RDM_files);


%% Format for RSA toolbox
clear RDMs
for i = 1:size(RDM_avg,3)
    RDMs(i).RDM = squeeze(RDM_avg(:,:,i));
    RDMs(i).name = clust_names{i};
    RDMs(i).color = rand(1,3);
end

%% Output T-values
magnitudes.mean_tvalues = means;
magnitudes.median_tvalues = medians;
magnitudes.dimensions = {'subjects';'categories';'ROI'};
magnitudes.categ_names = stim_names;
magnitudes.ROI_names = clust_names;


if ~dispFigs
    return
end


%% Distances
if exist('distances','var')
    figure('name','Distances','NumberTitle','off')
    bar(distances)
    xticklabel_rotate((1:numel(clust_names))-0.3,90,clust_names)
    % set(gca,'XTick',1:numel(clust_names))
    % set(gca,'XTickLabel',clust_names)
    ylabel('Distance (mm)')
    ylim([0 max(distances(:))+1])
    grid on
end

%% Mean RDM
plot_RDM(RDM_avg,'Mean',clust_names,stim_names,colmap)







%% misc
m = floor(sqrt(nclust));
n = ceil((nclust)/floor(sqrt(nclust)));
if ~exist('def_colors','var')
    if n_stims <= 7
        def_colors = lines(n_stims);
    else
        def_colors = hsv(n_stims);
    end
end


% %% T-values median
% figure('name','Median T-values','NumberTitle','off')
% miny = inf;
% maxy = -inf;
% for clust = 1:nclust
%     miny = min([miny;floor(min(min(mean(squeeze(medians(:,:,clust)),1) - max(std(squeeze(medians(:,:,clust)),[],1) ./ numel(RDM_files)))))]);
%     maxy = max([maxy;ceil(max(max(mean(squeeze(medians(:,:,clust)),1) + max(std(squeeze(medians(:,:,clust)),[],1)  ./ numel(RDM_files)))))]);
% end
% if miny > 0;miny = 0;end

% for clust = 1:nclust
%     subplot(m,n,clust)
%     plot_betas(mean(squeeze(medians(:,:,clust)),1),...
%         std(squeeze(medians(:,:,clust)),[],1) ./ numel(RDM_files),...
%         std(squeeze(medians(:,:,clust)),[],1) ./ numel(RDM_files),...
%         'Mean',stim_names,def_colors)

%     title(clust_names{clust})
%     ylim([miny maxy])
%     ylabel('Median T-value')
% end



% %% T-values median, mean of all clusters
% figure('name','Median T-values (all clusters)','NumberTitle','off')
% plot_betas(mean(mean(medians,3),1),...
%     std(mean(medians,3),[],1) ./ numel(RDM_files),...
%     std(mean(medians,3),[],1) ./ numel(RDM_files),...
%     'Mean',stim_names,def_colors)

% title('Mean of all clusters')
% ylim([miny maxy])
% ylabel('Median T-value')





%% T-values mean
figure('name','Mean T-values','NumberTitle','off')
miny = inf;
maxy = -inf;
for clust = 1:nclust
    miny = min([miny;floor(min(min(mean(squeeze(means(:,:,clust)),1) - max(std(squeeze(means(:,:,clust)),[],1) ./ numel(RDM_files)))))]);
    maxy = max([maxy;ceil(max(max(mean(squeeze(means(:,:,clust)),1) + max(std(squeeze(means(:,:,clust)),[],1) ./ numel(RDM_files)))))]);
end
if miny > 0;miny = 0;end
for clust = 1:nclust
    subplot(m,n,clust)
    plot_betas(mean(squeeze(means(:,:,clust)),1),...
        std(squeeze(means(:,:,clust)),[],1) ./ numel(RDM_files),...
        std(squeeze(means(:,:,clust)),[],1) ./ numel(RDM_files),...
        'Mean',stim_names,def_colors)

    title(clust_names{clust})
    ylim([miny maxy])
    ylabel('Mean T-value')
end




%% T-values mean, mean of all clusters
figure('name','Mean T-values (all clusters)','NumberTitle','off')
plot_betas(mean(mean(means,3),1),...
    std(mean(means,3),[],1) ./ numel(RDM_files),...
    std(mean(means,3),[],1) ./ numel(RDM_files),...
    'Mean',stim_names,def_colors)

title('Mean of all clusters')
ylim([miny maxy])
ylabel('Mean T-value')


%% Correlations
if nclust > 1
    RDMcorr(RDM_avg,nclust);
end






%% Mean of all clusters
plot_RDM(mean(RDM_avg,3),'Mean of all clusters',{'Mean of all clusters'},stim_names,colmap)













