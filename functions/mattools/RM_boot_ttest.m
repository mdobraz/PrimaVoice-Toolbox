function [pvalue, bootmean] = RM_boot_ttest(X,Y,iterations,tail)
% Repeated measures permutation test
% [pvalue, bootmean]  = boot_ttest(X,Y,iterations,tail), performs a
% bootstrapped test of the null hypothesis H0 that the data in vectors
% X and Y come from the same distribution.
%
% OUTPUT is the pvalue, the probability, under the null hypothesis, of
% observing a value as extreme or more extreme of the test statistic.
%
% iterations is the number of permutations, default: 1000
%
% There are three options for tail: 
% * 'both' — Means are not equal (two-tailed test). This is the default, when tail is unspecified.
% * 'right' — Mean of X is greater than mean of Y (right-tail test)
% * 'left' — Mean of X is less than mean of Y (left-tail test)

%% input arguments
if (nargin < 2)
  error(' There must be at least two args, sampleref and sampletest ');
elseif (nargin == 2)
   iterations = 1000;
   tail = 'both';
elseif (nargin == 3)
    tail = 'both';
elseif (nargin > 4)
    error(' Too many arguments !!!! ');
end


%% Transpose if necessary
if size(X,2) == 1; X = X';end
if size(Y,2) == 1; Y = Y';end

%% Check size of X,Y
if length(X) ~= length(Y)
    error('Vectors must be the same length')
end

%% Bootstrapping
bootmean = zeros(iterations,1);
Z = [X;Y];
for i = 1:iterations
    perm = logical(randi(2,1,length(X)) - 1);
    Zperm = [flipud(Z(:,perm)) Z(:,~perm)];
    bootmean(i) = diff(mean(Zperm,2));
end

%% pvalue
meandiff = mean(X) - mean(Y);
pcount = 0;
for it = 1:iterations
    if strcmp(tail,'both')
        if abs(bootmean(it)) >= abs(meandiff)
            pcount = pcount + 1;
        end
    elseif strcmp(tail,'left')
        if bootmean(it) <= meandiff
            pcount = pcount + 1;
        end 
    elseif strcmp(tail,'right')
        if bootmean(it) >= meandiff
            pcount = pcount + 1;
        end 
    end
end
pvalue = (pcount + 1) / (iterations + 1);


%% Diff of means
function DM = diff_mean(x,index)

y1 = x(1:index);
y2 = x(index+1:end);

DM = mean(y1) - mean(y2);