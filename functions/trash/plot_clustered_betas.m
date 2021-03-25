%% Options
% p value
p_val = p_val_peak(1);

% peak (0) or cluster (1) threshold
extent_thres = 0;

% vs_silence betas or not (comment if not)
name_suffix = '_vs-silence';


max_n_clust = 8;
confidence_level = 0.95;
nconds = contrasts.nconds;

conds = 2:nconds; % vs_silence conds
% conds = nconds+1:nconds+nconds-1; % betas conds



% % all contrasts of interest
% pref_cons = nconds+nconds+1:length(contrasts.names); % will determine preferences based on the last contrasts
% pref_cons  = [1 pref_cons ]; % add sound_vs_silence contrast

% sound_vs_silence only
pref_cons = 17;



def_colors = get(groot,'defaultAxesColorOrder');
line_space = 0.05;
line_col = [0.5 0.5 0.5];





%%%%% ELOUK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conds = [2 3 5];
% def_colors(3,:) = def_colors(4,:);




%% t-tests conditions
cond_triang = tril(~diag(conds>0));
first_diag = ~tril(cond_triang,-2) & cond_triang;
first_diag_i = nan(sum(first_diag(:)),2);
[first_diag_i(:,1),first_diag_i(:,2)] = ind2sub(size(first_diag),find(first_diag));

mid_conds = cond_triang & ~ first_diag & ~tril(cond_triang,-size(cond_triang,1)+1);
mid_conds_i = nan(sum(mid_conds(:)),2);
[mid_conds_i(:,1),mid_conds_i(:,2)] = ind2sub(size(mid_conds),find(mid_conds));

last_cond = tril(cond_triang,-size(cond_triang,1)+1);
last_cond_i = nan(sum(last_cond(:)),2);
[last_cond_i(:,1),last_cond_i(:,2)] = ind2sub(size(last_cond),find(last_cond));















% cond_triang = tril(~diag(conds>0));
% first_diag = ~tril(cond_triang,-2) & cond_triang;
% first_diag_i = nan(sum(first_diag(:)),2);
% [first_diag_i(:,1),first_diag_i(:,2)] = ind2sub(size(first_diag),find(first_diag));
% if length(conds) > 3
%     mid_conds = cond_triang & ~ first_diag & ~tril(cond_triang,-size(cond_triang,1)+1);
%     mid_conds_i = nan(sum(mid_conds(:)),2);
%     [mid_conds_i(:,1),mid_conds_i(:,2)] = ind2sub(size(mid_conds),find(mid_conds));
%     last_cond = tril(cond_triang,-size(cond_triang,1)+1);
%     last_cond_i = nan(sum(last_cond(:)),2);
%     [last_cond_i(:,1),last_cond_i(:,2)] = ind2sub(size(last_cond),find(last_cond));
% elseif length(conds) == 3
%     mid_conds_i = zeros(0,2);
%     last_cond_i = [3 1];
% else
%     mid_conds_i = zeros(0,2);
%     last_cond_i = zeros(0,2);
% end
% 
% ttest_conds = [first_diag_i;mid_conds_i;last_cond_i];











% cond_triang = tril(~diag(conds>0));
% ttest_conds = nan(sum(cond_triang(:)),2);
% ttest_conds_level = nan(sum(cond_triang(:)),1);
% 
% i1 = 1;
% for d = 1:length(conds)-1
%     this_diag = diag(conds(1:end-d)>0,-d);
%     this_diag_i = nan(sum(this_diag(:)),2);
%     [this_diag_i(:,1),this_diag_i(:,2)] = ind2sub(size(this_diag),find(this_diag));
%     i2 = i1 + size(this_diag_i,1) - 1;
%     ttest_conds(i1:i2,:) = this_diag_i;
%     ttest_conds_level(i1:i2) = d;
%     i1 = i2 + 1;
% end





p_val_str = sprintf('%g',p_val);
p_val_str = strrep(p_val_str,'.','p');

if extent_thres
    thres_suffix = [p_val_str '_clust-thres'];
else
    thres_suffix = [p_val_str '_peak-thres'];
end

if ~exist('name_suffix','var')
    name_suffix = thres_suffix;
else
    name_suffix = [thres_suffix name_suffix];
end




for icon = 1:length(pref_cons)
    con = pref_cons(icon);
    fprintf('%s\n',contrasts.names{con})
    clust_betas_file = fullfile(paths.results_multi,['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_clusters_betas_' name_suffix '.mat']);
    if exist(clust_betas_file,'file')
        load(clust_betas_file)
        nclust = length(cluster.maxT);
        if nclust > max_n_clust
            nclust = max_n_clust;
        end
        figure('name',contrasts.names{con},'NumberTitle','off')
        m = floor(sqrt(nclust));
        n = ceil((nclust)/floor(sqrt(nclust)));
        for clust = 1:nclust
            L = nan(length(conds),1);
            med = nan(length(conds),1);
            U = nan(length(conds),1);
%             PX = nan(length(conds),1);
            
            subplot(m,n,clust)
            hold on
            for i = 1:length(conds)
                [~,~,med(i),L(i),U(i)] = CI_values(cluster.xpd{clust,i},cluster.csypd{clust,i},confidence_level);
            end
            set(gca,'XTick',conds)
            ticks = cell(length(conds),1);
            for i = 1:length(conds)
                C = strsplit(contrasts.names{conds(i)},'_');
                ticks{i} = C{1};
            end
            set(gca,'XTickLabel',ticks)
            xlabel('Condition')
            ylabel('Activation')
            title(sprintf('Clust %i, nvox: %i, max T = %1.2f',clust,numel(cluster.values{clust,i}),cluster.maxT(clust)))
%             coords = cluster.max_coord{clust};
%             title(sprintf('Clust %i, nvox: %i, max T = %1.2f at %i, %i, %i',clust,numel(cluster.values{clust,i}),cluster.maxT(clust),coords(1),coords(2),coords(3)))
%             bar(conds,med,'FaceColor',def_colors(1,:),'EdgeColor','k');
            
            for b = 1:length(conds)
                bar(conds(b),med(b),'FaceColor',def_colors(b,:),'EdgeColor','k');
            end
            errorbar(conds,med,L,U,'.k','LineWidth',2)
            
            top = max(med + U);
            bot = min(med - L);
            st = (top - bot) / 30;
            
            
%             for i = 1:size(ttest_conds,1)
%                 c1 = ttest_conds(i,1);
%                 c2 = ttest_conds(i,2);
%                 [pvalue,ptext] = boot_ttest(cluster.values{clust,c1},cluster.values{clust,c2},1000,'both','median');
%                 if pvalue < 0.05
%                     yline = top + (st * ttest_conds_level(i));
%                     x1line = (min([conds(c1) conds(c2)]) + line_space);
%                     x2line = (max([conds(c1) conds(c2)]) - line_space);
%                     plot([x1line x2line],[yline yline],'color',line_col)
%                     ytext = yline + st/5;
%                     text(mean([conds(c1) conds(c2)]),ytext,ptext,'color',line_col,'HorizontalAlignment','center')
%                 end
%             end
            
            
            
            for i = 1:size(first_diag_i,1)
                c1 = first_diag_i(i,1);
                c2 = first_diag_i(i,2);
                [pvalue,ptext] = boot_ttest(cluster.values{clust,c1},cluster.values{clust,c2},1000,'both','mean');
                if pvalue < 0.05
                    yline = top + st;
                    x1line = (min([conds(c1) conds(c2)]) + line_space);
                    x2line = (max([conds(c1) conds(c2)]) - line_space);
                    plot([x1line x2line],[yline yline],'color',line_col)
                    ytext = yline + st/5;
                    text(mean([conds(c1) conds(c2)]),ytext,ptext,'color',line_col,'HorizontalAlignment','center')
                end
            end
            
            shift = 0;
            for i = 1:size(mid_conds_i,1)
                c1 = mid_conds_i(i,1);
                c2 = mid_conds_i(i,2);
                [pvalue,ptext] = boot_ttest(cluster.values{clust,c1},cluster.values{clust,c2},1000,'both','mean');
                if pvalue < 0.05
                    yline = top + (abs(c1-c2) * st) + (shift * st);
                    x1line = (min([conds(c1) conds(c2)]) + line_space);
                    x2line = (max([conds(c1) conds(c2)]) - line_space);
                    plot([x1line x2line],[yline yline],'color',line_col)
                    ytext = yline + st/5;
                    text(mean([conds(c1) conds(c2)]),ytext,ptext,'color',line_col,'HorizontalAlignment','center')
                    shift = shift + 1;
                end
            end
            
            c1 = last_cond_i(1,1);
            c2 = last_cond_i(1,2);
            [pvalue,ptext] = boot_ttest(cluster.values{clust,c1},cluster.values{clust,c2},1000,'both','mean');
            if pvalue < 0.05
                yline = top + (abs(c1-c2) * st) + ((shift-1) * st);
                x1line = (min([conds(c1) conds(c2)]) + line_space);
                x2line = (max([conds(c1) conds(c2)]) - line_space);
                plot([x1line x2line],[yline yline],'color',line_col)
                ytext = yline + st/5;
                text(mean([conds(c1) conds(c2)]),ytext,ptext,'color',line_col,'HorizontalAlignment','center')
            end
            
            
        end
    end
end




