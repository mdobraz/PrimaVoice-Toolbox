function NormLoudness(stims_dir,FS_desired,time_shifting)
% Normalize the loudness stimuli contained in the 'stims_dir' folder
% Uses the Zwicker and Fastl (1999) model for loudness for time-varying sounds
% Loudness calculation is based both on the max loudness and a longer term
% loudness (see LT_ratio, by default the loudest 10% of the signal)
% The loudness estimate is a weighted average of those two values, weighted
% by the parameter 'max_to_LT_ratio', 0.75 by default, 1 meaning only the
% max is used, 0 meaning only the long term loudness is used.
%
% The normalization is done iteratively until the mean loudness ratio -1
% between all stimuli is < 'precision', 0.001 by default
%
% Normalized stimuli are save in a new folder 'stims_dir_Loudness_Normalized'
%
% If 'FS_desired' is specified, the stimuli will also be resampled to the
% desired FS

%% editable variables
stim_extension = 'wav';
LT_ratio = 20; % long term ratio in %. Percentage of the signal used to compute Nx (Nx: loudness exceeded during x percent of the signal)
max_to_LT_ratio = 0.65; % weight of the max vs LT ratio (0 to 1)
precision = 0.005;
new_stims_dir = sprintf('%s_Loudness_Normalized',stims_dir);


%% Copy files to new folder & get stimuli filenames
if exist(new_stims_dir,'dir')
    delete(sprintf('%s/*',new_stims_dir)) % delete content of folder
else
    mkdir(new_stims_dir) % create folder if non-existant
end

copyfile(stims_dir,new_stims_dir);

stimuli = dir(sprintf('%s/*.%s',stims_dir,stim_extension));
n_stims = length(stimuli);


% stims = dir(sprintf('%s/*.wav',stims_dir));
% n_stims = 0;
% for i = 1:length(stims) % exclude '.' and '..'
%     if length(stims(i).name) > 4
%         if strcmp(stims(i).name(end-2:end),stim_extension)
%             n_stims = n_stims + 1;
%             stimuli(n_stims,:) = stims(i).name;
%         end
%     end
% end


%% first Loudness calculation
for i = 1:n_stims
    [signal,FS] = audioread(sprintf('%s/%s',new_stims_dir,stimuli(i).name)); % get signal
    L = Loudness_TimeVaryingSound_Zwicker(signal,FS,'mic','free',LT_ratio);
    eval(sprintf('L%i = L;',i));
end

% max
maxL = nan(n_stims,1);
for i = 1:n_stims
    eval(sprintf('maxL(i) = max(L%i.InstantaneousLoudness);',i,i));
end

% Nx
NxL = nan(n_stims,1);
for i = 1:n_stims
    eval(sprintf('NxL(i) = L%i.Nx;',i,i));
end

% mean
Loudness = (maxL .* max_to_LT_ratio) + (NxL .* (1 - max_to_LT_ratio)); % weighted mean of max and Nx
mean_Loudness = mean(Loudness);
% mean_Loudness = min(Loudness);

fact = mean_Loudness ./ Loudness;

%% Iterations
iteration = 0;
% fact = zeros(n_stims,1);
while mean(abs(fact-1)) > precision
    iteration = iteration + 1;
    fprintf('Iteration %i, loudness ratio = %1.3f (precision = %1.3f)\n',iteration,mean(abs(fact-1)),precision)
    for i = 1:n_stims
        filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
        [signal,FS] = audioread(filename); % get signal
        signal = signal .* fact(i);
        L = Loudness_TimeVaryingSound_Zwicker(signal,FS,'mic','free',LT_ratio);
        eval(sprintf('new_L%i = L;',i));
        audiowrite(filename,signal,FS)
    end
    
    % max
    maxL = nan(n_stims,1);
    for i = 1:n_stims
        eval(sprintf('maxL(i) = max(new_L%i.InstantaneousLoudness);',i,i));
    end

    % Nx
    NxL = nan(n_stims,1);
    for i = 1:n_stims
        eval(sprintf('NxL(i) = new_L%i.Nx;',i,i));
    end

    % mean
    Loudness = (maxL .* max_to_LT_ratio) + (NxL .* (1 - max_to_LT_ratio)); % weighted mean of max and Nx
%     mean_Loudness = mean(Loudness);
    mean_Loudness = min(Loudness);

    fact = mean_Loudness ./ Loudness;
end

fprintf('Last iteration, loudness ratio = %1.3f (precision = %1.3f)\n',mean(abs(fact-1)),precision)

%% Time shifting
if exist('time_shifting','var')
    if time_shifting
        for i = 1:n_stims
            eval(sprintf('new_Sones(%i,1:length(new_L%i.InstantaneousLoudness)) = new_L%i.InstantaneousLoudness;',i,i,i));
        end

        lag = nan(n_stims,1);
        for i = 1:n_stims
            [c,lags] = xcorr(new_Sones(i,:),mean(new_Sones));
            [~ ,b] = max(c);
            lag(i) = lags(b);
        end

        lag = max(lag) - lag;

        fprintf('Time adjustments')
        length_signal = nan(n_stims,1);
        for i = 1:n_stims % add zeros at the beginning
            time_lag = L.time(lag(i)+1);
            filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
            [signal,FS] = audioread(filename); % get signal
            sample_lag = round(time_lag * FS);
            if sample_lag
                signal = vertcat(zeros(sample_lag,1),signal);
                audiowrite(filename,signal,FS)
            end
            length_signal(i) = length(signal);
        end

        [max_length_signal,longest_stim] = max(length_signal);
        for i = 1:n_stims % add zeros at the end
            fprintf('.')

            filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
            [signal,FS] = audioread(filename); % get signal
            if i ~= longest_stim
                sample_lag = max_length_signal - length_signal(i);
                signal = vertcat(signal,zeros(sample_lag,1));
            end
            L = Loudness_TimeVaryingSound_Zwicker(signal,FS,'mic','free',LT_ratio);
            eval(sprintf('new_L%i = L;',i));
            audiowrite(filename,signal,FS)
        end
        fprintf('\nSignal length = %i samples\n',max_length_signal)
    end
end


%% Plot
for i = 1:n_stims
    eval(sprintf('Levels(%i,1:length(L%i.InstantaneousLoudnessLevel)) = L%i.InstantaneousLoudnessLevel;',i,i,i));
    eval(sprintf('new_Levels(%i,1:length(new_L%i.InstantaneousLoudnessLevel)) = new_L%i.InstantaneousLoudnessLevel;',i,i,i));
    eval(sprintf('Sones(%i,1:length(L%i.InstantaneousLoudness)) = L%i.InstantaneousLoudness;',i,i,i));
    eval(sprintf('new_Sones(%i,1:length(new_L%i.InstantaneousLoudness)) = new_L%i.InstantaneousLoudness;',i,i,i));
end

figure
subplot(2,2,1)
plot(Levels')
title('Before')
axis([0 length(Levels) 0 max([max(Levels(:)) max(new_Levels(:))])+5])
ylabel('Phones')
subplot(2,2,2)
plot(new_Levels')
title('After')
axis([0 length(Levels) 0 max([max(Levels(:)) max(new_Levels(:))])+5])
subplot(2,2,3)
plot(Sones')
% title('before')
axis([0 length(Sones) 0 max([max(Sones(:)) max(new_Sones(:))])+2])
ylabel('Sones')
xlabel('Time (ms)')
subplot(2,2,4)
plot(new_Sones')
axis([0 length(Sones) 0 max([max(Sones(:)) max(new_Sones(:))])+2])
xlabel('Time (ms)')
% title('after')
hold off


%% Resampling
if exist('FS_desired','var')
    fprintf('Resampled to %i Hz\n',round(FS_desired))
    for i = 1:n_stims
        filename = sprintf('%s/%s',new_stims_dir,stimuli(i).name);
        [signal,FS] = audioread(filename); % get signal
        signal = resample(signal,round(FS_desired),round(FS)); % resample the signal to the desired sample rate 
        audiowrite(filename,signal,round(FS_desired))
    end
end



%% Sound original and new mixtures
% clear mixture
% for i = 1:n_stims
%     [signal,FS] = audioread(sprintf('%s/%s',stims_dir,stimuli(i,:))); % get signal
%     if ~exist('mixture','var')
%         mixture = signal;
%     else
%         mixture = mixture + signal;
%     end
% end
% 
% clear new_mixture
% for i = 1:n_stims
%     [signal,FS] = audioread(sprintf('%s/%s',new_stims_dir,stimuli(i,:))); % get signal
%     if ~exist('new_mixture','var')
%         new_mixture = signal;
%     else
%         new_mixture = new_mixture + signal;
%     end
% end
% 
% sound(mixture,FS)
% sound(new_mixture,FS)
