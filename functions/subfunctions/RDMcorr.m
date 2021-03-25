function corrMat = RDMcorr(RDM,nclust)

cond_triang = tril(~diag(1:nclust));
pairs = nan(sum(cond_triang(:)),2);
[pairs(:,1),pairs(:,2)] = ind2sub(size(cond_triang),find(cond_triang));

cond_triang = tril(~diag(1:size(RDM,1)));
corrMat = zeros(nclust,nclust);

for p = 1:size(pairs,1)
    c1 = squeeze(RDM(:,:,pairs(p,1)));
    c1 = c1(cond_triang(:));

    c2 = squeeze(RDM(:,:,pairs(p,2)));
    c2 = c2(cond_triang(:));

    corrMat(pairs(p,1),pairs(p,2)) = corr(c1,c2,'type','Spearman');
end

corrMat = corrMat + corrMat' + eye(nclust,nclust);

figure('name','Spearman correlation between clusters','NumberTitle','off')

imagesc(corrMat,[-1 1]);
import rsa.fig.*
cols=colorScale([0 0.5 1; 0.5 0.5 0.5; 1 0 0],256);
colormap(cols);
colorbar
axis square
title('Spearman correlation between clusters')
set(gca,'XTick',1:nclust)
set(gca,'YTick',1:nclust)
set(gca,'TickLength',[0 0])









