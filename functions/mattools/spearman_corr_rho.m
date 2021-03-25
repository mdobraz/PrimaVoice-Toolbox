function rho = spearman_corr_rho(X,index)

if size(X,1) == 1 || size(X,2) == 1 % if vector

	y1 = X(1:index);
	y2 = X(index+1:end);

	% Transpose if necessary
	if size(y1,1) == 1; y1 = y1';end
	if size(y2,1) == 1; y2 = y2';end

	rho = corr([y1 y2],'type','Spearman');
elseif size(X,2) == 2 % if n-by-2 matrix
	rho = corr(X,'type','Spearman');
else
	error('X must either be a vector or a n-by-2 matrix')
end

rho = rho(1,2);
