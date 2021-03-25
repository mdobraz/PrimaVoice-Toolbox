function plot_RDM(RDM,figname,RDM_names,stim_labels,colmap)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if strcmp(colmap,'rsa')
    colmap = RDMcolormap;
end


nclust = numel(RDM_names);
n_stims = size(RDM,1);

if nclust > 1
    m = floor(sqrt(nclust+1));
    n = ceil((nclust+1)/floor(sqrt(nclust+1)));
end

%% RDMs
figure('name',[figname ' brain RDM'],'NumberTitle','off')
for clust = 1:nclust
    if nclust > 1
        subplot(m,n,clust)
    end
    % imagesc(squeeze(RDM(:,:,clust)))

    thisRDM = squeeze(RDM(:,:,clust));
    alpha = ~isnan(thisRDM);
    image(scale01(rankTransform_equalsStayEqual(thisRDM,1)),'CDataMapping','scaled','AlphaData',alpha);
    set(gca,'CLim',[0 1],'CLimMode','manual');
    colormap(gca,colmap);    
    title(RDM_names{clust})
    axis square %off;


    if exist('stim_labels','var') && ~isempty(stim_labels)
        if iscell(stim_labels) || isvector(stim_labels)
            set(gca,'XTick',1:length(stim_labels))
            
            if iscell(stim_labels)
                if length(stim_labels{1}) < 3
                    set(gca,'XTickLabel',deunderscore(stim_labels))
                elseif length(stim_labels{1}) < 8
                    xticklabel_rotate(1:length(stim_labels),90,deunderscore(stim_labels))
                end
            else
                set(gca,'XTickLabel',stim_labels)
            end
            set(gca,'YTick',1:length(stim_labels))
            set(gca,'YTickLabel',deunderscore(stim_labels))
            set(gca,'TickLength',[0 0])
        elseif isstruct(stim_labels)
            addImageSequenceToAxes(gca,stim_labels); 
        end
    else
        set(gca,'XTick',[],'YTick',[]);
    end
end

if nclust > 1
    subplot(m,n,nclust+1); cla;
    imagesc(rankTransform(squareRDM(thisRDM),1),[0 1]);  cla;
    % ht=text(n/2,n/2,{['\bfeach dissimilarity matrix (',num2str(n_stims),'^2)'], 'separately rank-transformed', 'and scaled into [0,1]'},'HorizontalAlignment','Center','FontUnits','normalized');
    % set(ht,'FontSize',.06);
    axis square off;
    colormap(gca,colmap);
end
colorbar

