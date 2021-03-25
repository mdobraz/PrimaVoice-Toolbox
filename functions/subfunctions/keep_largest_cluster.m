function X = keep_largest_cluster(X,conn,prob)
% Keep largest cluster from an image

if exist('prob','var')
    [labels,nclust] = bwlabeln(X > prob,conn); % 6-connected, 18-connected or 26-connected
else
    [labels,nclust] = bwlabeln(X,conn); % 6-connected, 18-connected or 26-connected
end


clust_sizes = nan(nclust,1);
for clust = 1:nclust
    clust_sizes(clust) = sum(labels(:) == clust);
end
[~,I] = max(clust_sizes);
X(labels(:) ~= I) = 0;



