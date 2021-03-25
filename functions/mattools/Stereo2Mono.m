function Stereo2Mono(stims_dir)
% stimuli will be lengthen or shorten to the desired NS (number of samples)

%% editable variables
stim_extension = 'wav';
new_stims_dir = sprintf('%s_Mono',stims_dir);


stimuli = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stimuli);

%% Resampling
%fprintf('Resampled to %i Hz\n',round(FS_desired))
j = 0;
I = [];
for i = 1:n_stims
    filename = sprintf('%s/%s',stims_dir,stimuli(i).name);
    ai = audioinfo(filename);
    if ai.NumChannels > 1
        j = j + 1;
        I(end+1) = i;
    end
end

if j > 0
    %% Copy files to new folder & get stimuli filenames
    if exist(new_stims_dir,'dir')
        delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
    else
        mkdir(new_stims_dir) % create folder if non-existant
    end

    copyfile(stims_dir,new_stims_dir);


    for i = I
        filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
        [signal,FS] = audioread(filename); % get signal
        signal = mean(signal,2);
        audiowrite(filename,signal,FS)
    end
    fprintf('\n%i files were converted to Mono\n\n',j)
else
    fprintf('\nAll files were already Mono, nothing was done\n\n')
end
    
    
    
    
    
    