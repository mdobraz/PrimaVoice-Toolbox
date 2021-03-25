function plot_betas(X,L,U,mask_contrast,stim_labels,def_colors)

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

hold on
for i = 1:numel(stim_labels)
    bar(i,X(i),'FaceColor',def_colors(i,:),'EdgeColor','k');
end
errorbar(1:numel(stim_labels),X,L,U,'.k','LineWidth',2)

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

hold off
% title(sprintf('Cluster %i',clust))
% ylim([miny maxy])
% ylabel('Median T-value')
