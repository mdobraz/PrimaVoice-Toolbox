function [pvalue,ptext,bootdist] = boot_ttest(X,Y,iterations,tail,which_stat)
% Regis March 2011
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
% * 'both' ? Means are not equal (two-tailed test). This is the default, when tail is unspecified.
% * 'right' ? Mean of X is greater than mean of Y (right-tail test)
% * 'left' ? Mean of X is less than mean of Y (left-tail test)

%% input arguments
if (nargin < 2)
  error(' There must be at least two args, sampleref and sampletest ');
elseif (nargin > 5)
    error(' Too many arguments !!!! ');
end

if ~exist('iterations','var')
    iterations = 1000;
end

if ~exist('tail','var')
    tail = 'both';
end

if ~exist('which_stat','var')
    which_stat = 'mean';
end


%% Transpose if necessary
if size(X,2) == 1; X = X';end
if size(Y,2) == 1; Y = Y';end

%% Bootstrapping
if strcmp(which_stat,'mean')
    bootdist = bootstrp(iterations,@diff_mean,[X Y],length(X));
    statdiff = mean(X) - mean(Y);
elseif strcmp(which_stat,'median')
    bootdist = bootstrp(iterations,@diff_median,[X Y],length(X));
    statdiff = median(X) - median(Y);
else
    error('which_stat must be ''mean'' or ''median''')
end



%% pvalue
if statdiff < min(bootdist) || statdiff > max(bootdist)
    if strcmp(tail,'both')
        pvalue = 0;
    elseif strcmp(tail,'left')
        if statdiff < min(bootdist)
            pvalue = 0;
        else
            pvalue = 1;
        end
    elseif strcmp(tail,'right')
        if statdiff > max(bootdist)
            pvalue = 0;
        else
            pvalue = 1;
        end
    end
    if pvalue == 0
        fprintf('p-value < %g\n',1/iterations)
    else
        fprintf('p-value = 1\n')
    end
else
    nbins = 1000;
    pd_D = fitdist(bootdist,'kernel');
    xpd = min(bootdist):1/nbins:max(bootdist);
    ypd = pdf(pd_D,xpd);
    csypd = cumsum(ypd);
    csypd = csypd ./ max(csypd);
    [~,yprob] = min(abs(xpd - statdiff));
    prob = csypd(yprob);
    
    if strcmp(tail,'both')
        if prob > 0.5
            prob = 1- prob;
        end
        pvalue = prob * 2;
    elseif strcmp(tail,'left')
        pvalue = prob; 
    elseif strcmp(tail,'right')
        pvalue = 1- prob; 
    else
        error('tail must be ''both'', ''left'' or ''right''')
    end
end

if pvalue < 0.001
    ptext = '***';
elseif pvalue < 0.01
    ptext = '**';
elseif pvalue < 0.05
    ptext = '*';
else
    ptext = 'n.s';
end



% %% Diff of means
% function DM = diff_mean(x,index)

% y1 = x(1:index);
% y2 = x(index+1:end);

% DM = mean(y1) - mean(y2);

% %% Diff of means
% function DM = diff_median(x,index)

% y1 = x(1:index);
% y2 = x(index+1:end);

% DM = median(y1) - median(y2);