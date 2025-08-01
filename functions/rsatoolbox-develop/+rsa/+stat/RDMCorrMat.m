function corrMat=RDMCorrMat(RDMs,figPlotSpec,type)
% USAGE
% corrMat=RDMCorrMat(RDMs[,figPlotSpec])
%
% FUNCTION
% returns and optionally displays the correlation matrix (spearman) of a set
% of RDMs (can be square or upper triangle form and wrapped or bare).
%
% 08-2011 CW: slight modifications for how nans are handled
%__________________________________________________________________________
% Copyright (C) 2011 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if ~exist('type','var'),type='Spearman'; end; % Spearman rank should be default?

RDMs_bareVecs=unwrapRDMs(vectorizeRDMs(RDMs));%reduceRDMsToValidConditionSet(RDMs))); % This is of size [1 utv nRDMs]

[one,nRDMParams,nRDMs]=size(RDMs_bareVecs);

RDMs_cols=permute(RDMs_bareVecs,[2 3 1]); % This is of size [utv nRDMs (1)]

% For each pair of RDMs, ignore missing data only for this pair of RDMs
% (unlike just using corr, which would ignore it if ANY RDM had missing
% data at this point).
%corrMat=corrcoef(RDMs_cols)
if isequal(type,'Kendall_taua')
    for RDMI1 = 1:nRDMs
        for RDMI2 = 1 : nRDMs
            corrMat(RDMI1,RDMI2)=rankCorr_Kendall_taua(RDMs_cols(:,RDMI1), RDMs_cols(:,RDMI2));
        end
    end
else
    for RDMI1 = 1:nRDMs
        for RDMI2 = 1 : nRDMs
            corrMat(RDMI1,RDMI2)=corr(RDMs_cols(:,RDMI1), RDMs_cols(:,RDMI2), 'type', type, 'rows', 'complete');
        end
    end
end
    
for RDMI1 = 1:nRDMs
	corrMat(RDMI1,RDMI1) = 1; % make the diagonal artificially one
end

if exist('figPlotSpec','var')
    selectPlot(figPlotSpec);
    imagesc(corrMat,[-1 1]);
    cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
    colormap(cols); colorbar;
    axis square %off;
    title(['\bfRDM correlation matrix (',type,')'],'FontUnits','normalized','FontSize',1/max(nRDMs,30)*2);

	if isstruct(RDMs)
		for RDMI=1:nRDMs
			% text(RDMI,RDMI,RDMs(RDMI).name,'HorizontalAlignment','center','FontWeight','bold','Color',RDMs(RDMI).color,'FontUnits','normalized','FontSize',1/max(nRDMs,30));
            RDM_labels{RDMI} = RDMs(RDMI).name;
		end

        set(gca,'XTick',1:nRDMs)
        if length(RDM_labels{1}) < 3
            set(gca,'XTickLabel',deunderscore(RDM_labels))
        else
            xticklabel_rotate(1:length(RDM_labels),90,deunderscore(RDM_labels))
        end
        set(gca,'YTick',1:nRDMs)
        set(gca,'YTickLabel',deunderscore(RDM_labels))
        set(gca,'TickLength',[0 0])
	end
end


%% Ranked corrs
% figure('name','Ranked correlations','NumberTitle','off')
% [m n]=paneling(nRDMs,2/3);
if isstruct(RDMs)
    for i=1:nRDMs
        % subplot(m,n,i)
        figure('name','Ranked correlations','NumberTitle','off')
        title(RDMs(i).name);
        [R,I] = sort(corrMat(i,:),'descend');
        hold on
        for j = 2:nRDMs
            bar(j-1,R(j),'FaceColor',RDMs(I(j)).color,'EdgeColor','k');
        end
        set(gca,'XTick',1:nRDMs-1)
        if length(RDM_labels{1}) < 5
            set(gca,'XTickLabel',deunderscore(RDM_labels(I(2:end))))
        else
            xticklabel_rotate(1:nRDMs-1,90,deunderscore(RDM_labels(I(2:end))))
        end
        ylim([-1 1])
        hold off
    end
end

end%function
