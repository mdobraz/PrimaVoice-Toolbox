function Norm_max_RMS(stims_dir,percent_loudest)
% Normalize wav files in 'stims_dir' folder.
% Normalizatio is done by comparing the loudest portion of the signal
% 'percent_loudest'
% the grand average of all stims is used as a reference


stim_extension = 'wav';
new_stims_dir = sprintf('%s_Normalized',stims_dir);

% list sounds
stimuli = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stimuli);



% get min duration of all stims
min_length = inf;
FSs = nan(n_stims,1);
for i = 1:n_stims
    f = sprintf('s%i',i);
    [signal.(f),FSs(i)] = audioread(sprintf('%s/%s',stims_dir,stimuli(i).name)); % get signal
    if length(signal.(f)) < min_length
        min_length = length(signal.(f));
    end
end

% check FS
if length(unique(FSs)) > 1
    error('WAV files have different sampling rates')
else
    FS = unique(FSs);
end

% mean signal of all stims before min duration
sum_sig = zeros(min_length,1);
for i = 1:n_stims
    f = sprintf('s%i',i);
    sum_sig = sum_sig + signal.(f)(1:min_length);
end
sum_sig = sum_sig ./ max(sum_sig);

% find loudest term RMS (grand_level) of mean signal to use it as a reference (contains all frequencies encountered in all stims)
sum_sig2 = sum_sig .^ 2;
pks = findpeaks(sum_sig2);
pks = sqrt(pks);
LT = round(length(pks) * (percent_loudest / 100));
spks = sort(pks,'descend');
spks = spks(1:LT);
grand_level = mean(spks);

% create nez folder for normalized stims
if exist(new_stims_dir,'dir')
    delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
else
    mkdir(new_stims_dir) % create folder if non-existant
end
copyfile(stims_dir,new_stims_dir);

% normalize all stims
max_s = 0;
for i = 1:n_stims
    f = sprintf('s%i',i);
    s = signal.(f);
    s2 = s .^ 2;
    pks = findpeaks(s2);
    pks = sqrt(pks);
    LT = round(length(pks) * (percent_loudest / 100));
    spks = sort(pks,'descend');
    spks = spks(1:LT);
    level = mean(spks);
    s = s .* (grand_level / level);
    if max(abs(s)) > max_s
        max_s = max(abs(s));
    end
    signal.(f) = s;
end

% normalize to max to prevent clipping
for i = 1:n_stims
    f = sprintf('s%i',i);
    s = signal.(f);
    s = s ./ max_s;
    filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
    audiowrite(filename,s,FS)
end

% Ramped & Save grand average
signal = sum_sig; % get signal
signal = signal ./ max_s;
rampdur = 10;
rampsamples = round(FS * rampdur/1000);
t = linspace(0,pi,rampsamples);
offramp = (cos(t) + 1) / 2;
onramp = fliplr(offramp);
signal(1:rampsamples) = signal(1:rampsamples) .* onramp';
signal(end-rampsamples+1:end) = signal(end-rampsamples+1:end) .* offramp';
[~,a] = fileparts(stims_dir);
filename = sprintf('%s_Normalization_Reference.wav',a);
audiowrite(filename,signal,FS)



