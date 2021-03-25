function AdjustLength(stims_dir,NS_desired)
% stimuli will be lengthen or shorten to the desired NS (number of samples)

%% editable variables
stim_extension = 'wav';
new_stims_dir = sprintf('%s_Adjusted',stims_dir);


stimuli = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stimuli);

%% Adjusting
%fprintf('Length adjusted to %i samples\n',NS_desired)
nss = nan(n_stims,1);
bis = nan(n_stims,1);
for i = 1:n_stims
    filename = sprintf('%s/%s',stims_dir,stimuli(i).name);
    ai = audioinfo(filename);
    nss(i) = ai.TotalSamples;
    bis(i) = ai.BitsPerSample;
end

fprintf('min NS = %i, max NS = %i\n\n',min(nss),max(nss))

r = input('Adjust length (y/n)? ','s');

if strcmp(r,'y')

    %% Copy files to new folder & get stimuli filenames
    if exist(new_stims_dir,'dir')
        delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
    else
        mkdir(new_stims_dir) % create folder if non-existant
    end

    copyfile(stims_dir,new_stims_dir);


    for i = 1:n_stims
        if nss(i) ~= NS_desired
            filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
            [signal,FS] = audioread(filename); % get signal
            if nss(i) < NS_desired
                signal(end+1:NS_desired) = 0;
            else
                signal = signal(1:NS_desired);
            end
            audiowrite(filename,signal,FS)
        end
    end
end
    
    
    
    
    
    