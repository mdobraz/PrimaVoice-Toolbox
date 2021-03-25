function DM = diff_mean(x,index)
% Diff of means
y1 = x(1:index);
y2 = x(index+1:end);

DM = mean(y1) - mean(y2);

% %% Diff of means
% function DM = diff_median(x,index)

% y1 = x(1:index);
% y2 = x(index+1:end);

% DM = median(y1) - median(y2);