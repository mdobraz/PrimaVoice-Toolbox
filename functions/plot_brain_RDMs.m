
RDM_var_and_filename

RDMs = load(RDM_file);

%% nclust
nclust = length(RDMs.cluster);
m = floor(sqrt(nclust));
n = ceil((nclust)/floor(sqrt(nclust)));

%% Stim names
stim_names = RDMs.stim_names;


n_stims = size(RDMs.RDM,1);
if n_stims <= 7
    def_colors = lines(n_stims);
    % stim_labels = stim_names;
else
    def_colors = hsv(n_stims);
    % stim_labels = def_colors;
end

%% RDMs
clust_names = cell(nclust,1);
for clust = 1:nclust
    clust_names{clust} = sprintf('Cluster %i',clust);
end
plot_RDM(RDMs.RDM,[mask_contrast ' ' pdistance],clust_names,stim_names,colmap)


% return


%% T-values min max
miny = inf;
maxy = -inf;
for clust = 1:nclust
    miny = min([miny;floor(min(min([RDMs.cluster(clust).stim(:).Tvalues])))]);
    maxy = max([maxy;ceil(max(max([RDMs.cluster(clust).stim(:).Tvalues])))]);
end
if miny > 0;miny = 0;end


%% T-values median
figure('name',[mask_contrast ' Median T-values'],'NumberTitle','off')
for clust = 1:nclust
    subplot(m,n,clust)
    plot_betas([RDMs.cluster(clust).stim(:).medianT],...
        [RDMs.cluster(clust).stim(:).median_L],...
        [RDMs.cluster(clust).stim(:).median_U],...
        mask_contrast,stim_names,def_colors)

    title(sprintf('Cluster %i',clust))
    ylim([miny maxy])
    ylabel('Median T-value')
end

%% T-values mean
figure('name',[mask_contrast ' Mean T-values'],'NumberTitle','off')
for clust = 1:nclust
    subplot(m,n,clust)
    plot_betas([RDMs.cluster(clust).stim(:).meanT],...
        [RDMs.cluster(clust).stim(:).mean_L],...
        [RDMs.cluster(clust).stim(:).mean_U],...
        mask_contrast,stim_names,def_colors)

    title(sprintf('Cluster %i',clust))
    ylim([miny maxy])
    ylabel('Mean T-value')
end



%% T-values data points
figure('name',[mask_contrast ' T-values'],'NumberTitle','off')
for clust = 1:nclust
    subplot(m,n,clust)
    hold on
    med = [RDMs.cluster(clust).stim(:).medianT];
    ci1 = [RDMs.cluster(clust).stim(:).median_ci1];
    ci2 = [RDMs.cluster(clust).stim(:).median_ci2];
    

    xrnd = (rand(1,length(RDMs.cluster(clust).stim(1).Tvalues)) .* 0.8) - 0.4;
    for i = 1:numel(stim_names)
        Tvalues = RDMs.cluster(clust).stim(i).Tvalues;
        plot(xrnd + i,Tvalues,'.','Color',def_colors(i,:));
        % plot(xrnd + i,Tvalues,'.');
        plot([i-0.4 i+0.4],[med(i) med(i)],'k','LineWidth',2)
        % plot([i-0.4 i+0.4],[ci1(i) ci1(i)],':k','LineWidth',2)
        % plot([i-0.4 i+0.4],[ci2(i) ci2(i)],':k','LineWidth',2)
    end
    % errorbar(1:numel(stim_names),med,L,U,'.k','LineWidth',2)
    title(sprintf('Cluster %i',clust))
    set(gca,'XTick',1:length(stim_names))
    set(gca,'XTickLabel',stim_names)
    
    ylim([miny maxy])
    ylabel('T-value')
    hold off
end




%% Correlations
if nclust > 1
    RDMcorr(RDMs.RDM,nclust);
end




% %% Correlations
% if nclust > 1
%     cond_triang = tril(~diag(1:nclust));
%     pairs = nan(sum(cond_triang(:)),2);
%     [pairs(:,1),pairs(:,2)] = ind2sub(size(cond_triang),find(cond_triang));

%     cond_triang = tril(~diag(1:n_stims));
%     corrMat = zeros(nclust,nclust);

%     for p = 1:size(pairs,1)
%         c1 = squeeze(RDMs.RDM(:,:,pairs(p,1)));
%         c1 = c1(cond_triang(:));

%         c2 = squeeze(RDMs.RDM(:,:,pairs(p,2)));
%         c2 = c2(cond_triang(:));

%         corrMat(pairs(p,1),pairs(p,2)) = corr(c1,c2,'type','Spearman');
%     end

%     corrMat = corrMat + corrMat' + eye(nclust,nclust);

%     figure('name',[mask_contrast ' brain RDM norm cluster correlation'],'NumberTitle','off')

%     imagesc(corrMat,[-1 1]);
%     cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
%     colormap(cols);
%     colorbar
%     axis square
%     title('Spearman correlation between clusters')
%     set(gca,'XTick',1:nclust)
%     set(gca,'YTick',1:nclust)
% end













