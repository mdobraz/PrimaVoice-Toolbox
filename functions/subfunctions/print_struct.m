function print_struct(S,fid,sub_struct)

if nargin < 2
    fid = 1;
end
if nargin < 3
    sub_struct = 0;
end

for el = 1:numel(S)
    jnd = fieldnames(S(el));
    for i = 1:numel(jnd)
        if sub_struct
            fprintf(fid,'\n');
            eval(['fprintf(fid,''' repmat('    ',1,sub_struct) ''');'])
%             fprintf(fid,'\n\t');
        end
        fprintf(fid,'%s:',jnd{i});
        val = S(el).(jnd{i});
        if ischar(val)
            fprintf(fid,' %s',val);
        elseif iscell(val)
            if numel(val) == 1
                fprintf(fid,' %s',cell2str(val));
            else
                fprintf(fid,'\n');
                for c = 1:numel(val)
                    fprintf(fid,'\t%s\n',val{c});
                end
            end
        elseif isstruct(val)
            print_struct(val,fid,sub_struct+1)
        else
            if size(val,1) == 1
                fprintf(fid,' %g',val);
            else
                for r = 1:size(val,1)
                    fprintf(fid,'\n\t');
                    fprintf(fid,' %g',val(r,:));
                end
            end
        end
        if ~sub_struct
            fprintf(fid,'\n');
        end
    end
%     fprintf(fid,'\n');
end
fprintf(fid,'\n');