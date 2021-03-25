function verif_file = compute_T2s(ME_files,session,run,paths,BIDS)

n_TE = numel(ME_files);

%% Get json
ME_json = cell(n_TE,1);
echo_times = nan(1,n_TE);
for i = 1:n_TE
    ME_json{i} = get_corresp_json(ME_files{i});
    ME_infos = spm_jsonread(ME_json{i});
    echo_times(i) = ME_infos.EchoTime;
end

%% Get session, run
T2s.session = session;
T2s.run = run;

%% Get date, time
ME_infos = spm_jsonread(ME_json{1});
ME_AT = ME_infos.AcquisitionTime;

% Find corresponding bold run (because AcquisitionDateTime is missing in the echo json)
bold_runs = spm_BIDS(BIDS,'data','sub',paths.subject,'ses',sprintf('%02.0f',T2s.session),'type','bold');
bold_json = get_corresp_json(bold_runs{1});
bold_infos = spm_jsonread(bold_json);
bold_AT = bold_infos.AcquisitionDateTime;
Tidx = strfind(bold_AT,'T');
bold_date = bold_AT(1:Tidx);
T2s.scan_time = datetime([bold_date ME_AT],'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS');



%% Average each echo
if ~exist(paths.ME,'dir');mkdir(paths.ME);end % create ME folder if non-existant
avg_ME_files = cell(n_TE,1);
for i = 1:n_TE
    avg_ME_files{i} = fullfile(paths.ME,sprintf('averaged_echo%i.nii.gz',i));
    system(sprintf('%sfslmaths %s -Tmean %s',paths.FSL_prefix,ME_files{i},avg_ME_files{i})); % average merged anat runs
end

%% Merged averaged echos
merge_in_files = sprintf('%s ',avg_ME_files{:});
merge_out_file = fullfile(paths.ME,'averaged_echos_merged.nii.gz');
system(sprintf('%sfslmerge -t %s %s',paths.FSL_prefix,merge_out_file,merge_in_files)); % merge anat runs

% %% Mask brain
avg_file = avg_ME_files{1};
% Setting paths & -p option for bash scripts
debias_path = which('T1xT2BiasFieldCorrection.sh');
BET_path = which('T1xT2BET.sh');
IterREGBET_path = which('IterREGBET.sh');
if isfield(paths,'FSL_prefix') && ~isempty(paths.FSL_prefix)
    popt = ['-p ' paths.FSL_prefix];
else
    popt ='';
end

% Debias
in_file = avg_file;
system(sprintf('bash %s -t1 %s -t2 %s -bet 2 -f 0.1 %s',debias_path,in_file,in_file,popt));
debiased_file = fullfile(paths.ME,'averaged_echo1_debiased.nii.gz');

% BETing debiased average
in_file = debiased_file;
system(sprintf('bash %s -t1 %s -t2 %s -n 2 -f 0.1 %s',BET_path,in_file,in_file,popt));
BETed_file = fullfile(paths.ME,'averaged_echo1_debiased_BET.nii.gz');

% Apply anat brain mask to debiased averages (IterREGBET)
inw = debiased_file;
inb = BETed_file;
refb = paths.anat.brain;
system(sprintf('bash %s -inw %s -inb %s -refb %s -dof 6 -n 3 %s',IterREGBET_path,inw,inb,refb,popt));
brain_file = fullfile(paths.ME,'averaged_echo1_debiased_IRbrain.nii.gz');


% %%%%%%%%%%%%%%%%%%
% system(sprintf('%sfslmaths %s -thr 3500 -bin %s',paths.FSL_prefix,avg_ME_files{1},brain_file)); % merge anat runs
% %%%%%%%%%%%%%%%%%%

Pbrain = spm_vol(brain_file);
Ybrain = logical(spm_read_vols(Pbrain));


%% Computation
x = echo_times .* 1000; % echo times in ms

P = spm_vol(merge_out_file);
Echos = spm_read_vols(P);

% monoexponential function
y = @(b,x) b(1).*exp(-b(2).*x);


opts = optimset('MaxFunEvals',10000,'MaxIter',10000);
T2s_map = nan(size(Ybrain));
% S0map = nan(size(Yroi));

fprintf('\nComputing T2s_map, session %i, run %i...\n',T2s.session,T2s.run)
tic
for i = 1:size(Echos,1)
    for j = 1:size(Echos,2)
        for k = 1:size(Echos,3)
            if Ybrain(i,j,k) % if voxel is in the ROI mask
                yx = squeeze(Echos(i,j,k,:))'; % all echos
                [~,I] = sort(yx,'descend');
                if prod(I==1:n_TE)
                    OLS = @(b) sum((y(b,x)-yx).^2);
                    B = fminsearch(OLS,rand(2,1),opts);
                    T2s_map(i,j,k) = 1/B(2);
                else
                    T2s_map(i,j,k) = nan;
                end
                if T2s_map(i,j,k) > 100
                    T2s_map(i,j,k) = nan;
                end
%                 S0map(i,j,k)=B(1);
            end
        end
    end
end
toc


% T2s_map = Echos(:,:,:,1);
% T2s_map(~Ybrain) = nan;

T2s.mean_T2 = nanmean(T2s_map(:)) ./ 1000; % mean T2 of the ROI from ms to sec
T2s.mean_R = 1 / T2s.mean_T2; % in s-1

T2s_file = fullfile(paths.ME,sprintf('sub-%s_ses-%02.0f_run-%02.0f_T2s.mat',paths.subject,T2s.session,T2s.run));
save(T2s_file,'-struct','T2s');


T2map_file = fullfile(paths.ME,sprintf('sub-%s_ses-%02.0f_run-%02.0f_T2s-map.nii',paths.subject,T2s.session,T2s.run));
PT2 = P(1);
PT2.fname = T2map_file;
spm_write_vol(PT2,T2s_map);


% verif file
[ME_path,ME_name,ext] = fileparts(ME_files{1});
if strcmp(ext,'.gz'); [~,ME_name] = fileparts(ME_name); end
verif_file = fullfile(paths.ME,[ME_name '.nii.gz']);
movefile(brain_file,verif_file)

delete(fullfile(paths.ME,'averaged_*'))



