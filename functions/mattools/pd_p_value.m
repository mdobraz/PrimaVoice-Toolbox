function [pvalue,ptext] = pd_p_value(testvalue,xpd,csypd,iterations,tail)

if ~exist('tail','var')
    tail = 'both';
end

if testvalue < min(xpd) || testvalue > max(xpd)
    if strcmp(tail,'both')
        pvalue = 0;
    elseif strcmp(tail,'left')
        if testvalue < min(xpd)
            pvalue = 0;
        else
            pvalue = 1;
        end
    elseif strcmp(tail,'right')
        if testvalue > max(xpd)
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
    
    [~,yprob] = min(abs(xpd - testvalue));
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