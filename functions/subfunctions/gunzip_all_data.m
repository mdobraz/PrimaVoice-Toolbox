function BIDS = gunzip_all_data(BIDS,paths,current_subject)
% Uncompresses all nii files of a protocol
% if current_subject is set to 1, compresses only the data of the subject
% contained in paths.subject

opt_str = '';
f = '';

if nargin > 2
    if current_subject
        opt_str = [opt_str ',''sub'',''' paths.subject ''''];
    end
end

eval(['f = spm_BIDS(BIDS,''data''' opt_str ');']);

c = ~cellfun(@isempty,regexp(f,'\.gz$','ignorecase')); % cells ending by '.gz'
if any(c) % if any
    fprintf('Uncompressing %i file(s)...\n',sum(c))
    I = 1:length(c); % indexes of c
    for i = I(c) % indexes of cells ending by '.gz'
        fprintf('%s\n',f{i})
        system(['gunzip ' f{i}]); % use system command because matlab command keeps the original by default
    end
    BIDS = spm_BIDS(paths.dataset);
    fprintf('Done.\n')
end