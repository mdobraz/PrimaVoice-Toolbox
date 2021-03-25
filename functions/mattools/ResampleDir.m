function ResampleDir(stims_dir,FS_desired,nbits_desired)
% stimuli will be resampled to the desired FS

if ~exist('nbits_desired','var')
    nbits_desired = 16;
end

%% editable variables
stim_extension = 'wav';
new_stims_dir = sprintf('%s_Resampled',stims_dir);


%% Copy files to new folder & get stimuli filenames
if exist(new_stims_dir,'dir')
    delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
else
    mkdir(new_stims_dir) % create folder if non-existant
end

copyfile(stims_dir,new_stims_dir);

stimuli = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stimuli);

%% Resampling
fprintf('Resampling to %i Hz\n',round(FS_desired))
for i = 1:n_stims
    filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
    [signal,FS] = audioread(filename); % get signal    
    if FS ~= round(FS_desired)

       signal = resample(signal,round(FS_desired),round(FS)); % resample the signal to the desired sample rate 
        if max(abs(signal)) > 1
            signal = signal ./ max(abs(signal));
        end
        audiowrite(filename,signal,round(FS_desired),'BitsPerSample',nbits_desired)
    end
end