function DM = diff_median(x,index)
% Diff of medians

y1 = x(1:index);
y2 = x(index+1:end);

DM = median(y1) - median(y2);