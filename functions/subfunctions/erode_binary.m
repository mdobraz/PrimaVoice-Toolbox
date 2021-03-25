function X = erode_binary(X,type)
% erodes a binary vector X
% type = 'both' will do a standard erosion
% type = 'before' will erode only before a group of ones
% type = 'after'  will erode only after  a group of ones

X = double(X);

if nargin < 2
    type = 'both';
end

for i = 1:length(X)
    if X(i)
        if strcmp(type,'both') || strcmp(type,'before')
            if (i > 1) && ~X(i-1)
                X(i) = 2;
            end
        end
        if strcmp(type,'both') || strcmp(type,'after')
            if (i < length(X)) && ~X(i+1)
                X(i) = 2;
            end
        end
    end
end

X(X==2) = 0;
X = logical(X);