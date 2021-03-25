function RampOnOffset(stims_dir,rampdur)
% Add a cos ramp of rampdur (ms) to the onset and offset of sounds

%% editable variables
stim_extension = 'wav';

new_stims_dir = sprintf('%s_ramped_%ims',stims_dir,round(rampdur));


%% Copy files to new folder & get stimuli filenames
if exist(new_stims_dir,'dir')
    delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
else
    mkdir(new_stims_dir) % create folder if non-existant
end

copyfile(stims_dir,new_stims_dir);

stims = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stims);


for i = 1:n_stims
    filename = sprintf('%s/%s',new_stims_dir,stims(i).name);
    [signal,FS] = audioread(filename); % get signal
    rampsamples = round(FS * rampdur/1000);
    t = linspace(0,pi,rampsamples);
    offramp = (cos(t) + 1) / 2;
    onramp = fliplr(offramp);
    signal(1:rampsamples) = signal(1:rampsamples) .* onramp';
    signal(end-rampsamples+1:end) = signal(end-rampsamples+1:end) .* offramp';
    audiowrite(filename,signal,FS)
end

