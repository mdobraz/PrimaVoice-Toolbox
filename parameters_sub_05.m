clearvars('-except','script_id','recipient','HN','run_index','max_rank','this_params_file'); % clear all but 'recipient' and 'HN' which are used in RunScript & 'run_index' used for jackknife procedure

%% PARAMETERS
% All the parameters of the analysis are set here
% Run this script before starting any analysis or using any function. Re-run it if you made changes.

%% Paths
paths.dataset = '/hpc/banco/Primavoice_Data_and_Analysis'; % Directory where the data are (general BIDS folder containing all subjects)
paths.subject = '05'; % the subject you want to analyze (eg. if the folder containing the subject data is called 'sub-moumou', omit the prefix 'sub-' and leave only 'moumou')
paths.task = 'sparse'; % the task you want to analyze (eg. if the datafiles include 'task-sparse', omit the prefix 'task-' and leave only 'sparse')
paths.realign_name = 'spm_realign'; % general name of the analysis (to differentiate analyses that have different parameters)
paths.results_name = '8WM_9CSF_0mvt'; % to distinguish different analyses done of the same realigned data
% paths.results_sub_folder = 'R2s38'; % Usefull when computing resulst for a limited number of runs (e.g. jackknife runs). Leave blank '' otherwise.

% Main anatomic file(s)
paths.main_anat_suffix = 'T1w'; % Suffix of the main structural scan (last part of the anat scan filename [eg. 'T1w'])
paths.main_anat_session = 1; % Session at which the main structural scan was acquired.
paths.T2_anat_suffix = 'T2w'; % Suffix of the T2 structural scan (comment line if no T2 available)
paths.T2_anat_session = 1; % Session at which the T2w structural scan was acquired.

% Secondary (daily) anatomic file(s).
% paths.close2func_anat_suffix = 'T1wReor'; % Suffix (last part of the anat scan filename [eg. 'T1w'])
% 								 % of a structural scan that is spatially close to the functional scans (used for func to anat coregistration).
% 								 % This can be the suffix of only one scan, or of scans acquired before each functional session.
% 								 % Set to '' or comment if unavailable.
% paths.close2func_anat_session = 1; % Specify the session of the scan, if only one scan is available
% 								 % Set to '' or comment if anat scans were acquired before each functional session.

% Segmentation template (for non-human primates). Comment for humans.
paths.templates = '/hpc/banco/Primavoice_Data_and_Analysis/templates'; % folder where templates are stored
paths.template_name = 'MNI152'; % Template that will be used during the segmentation process. Name should match a 'case' in the switch found in functions/templates.m
                             % Edit file functions/templates.m to add a new template
% FSL & ANTS paths
paths.FSL_prefix = 'fsl5.0-'; % prefix to call FSL tools by the system (eg. on frioul, to call 'flirt' you have to type 'fsl5.0-flirt', so the prefix is 'fsl5.0-'). Leave '' if no prefix is needed
paths.ANTS_path = '/hpc/soft/ANTS/antsbin/bin'; % path to ANTS binaries
paths.ANTS_scripts = '/hpc/soft/ANTS/ANTs/Scripts'; % path to ANTS scripts
paths.brainvisa = '/hpc/soft/brainvisa/brainvisa/bin'; % path to brainvisa binaries

%% Session and runs
sessions = 'all'; % either a vector or 'all' if you want to analyze all sessions
% sessions = [2:14]; % either a vector or 'all' if you want to analyze all sessions
% sessions = [5 12]; % either a vector or 'all' if you want to analyze all sessions
% 
runs = 'all'; % either a vector (if you want the same run numbers in each session) or 'all' if you want to analyze every run of each session
% runs = 1:2; % either a vector (if you want the same run numbers in each session) or 'all' if you want to analyze every run of each session

% runs = 1:2;

% runs = [{6};...
%         {3}];

% to select particular runs in each session, do something like this (one cell per session, a vector in each cell):
% runs = [{1:3};...
%         {1:4}];


%%% To run the analysis on runs previously selected by R2s, motion or jackknife methods:
% filename = fullfile(fullfile(paths.dataset,['analysis_sub-' paths.subject]),'R2s_Run_selection_64.mat');
% % filename = fullfile(fullfile(paths.dataset,['analysis_sub-' paths.subject]),paths.realign_name,'MotionCorr_Run_selection_25.mat');
% % filename = fullfile(fullfile(paths.dataset,['analysis_sub-' paths.subject]),paths.realign_name,'MCxR2s_Run_selection_26.mat');
% % filename = fullfile(fullfile(paths.dataset,['analysis_sub-' paths.subject]),paths.realign_name,'Jackknife_Run_selection_32.mat');
% % filename = fullfile(fullfile(paths.dataset,['analysis_sub-' paths.subject]),paths.realign_name,'Jackknife_Run_selection_32_2_folds_group2.mat');
% load(filename,'selected_sessions','selected_runs')
% sessions = selected_sessions;
% runs = selected_runs;



%% Brain segmentation parameters
brain_seg.BETf = 0.5; % -f option of BET. (0->1); default=0.5; smaller values give larger brain outline estimates. 0.5 for human, 0.3 for macaques
brain_seg.flirt_cost = 'normmi'; % cost function for subject-template registration with flirt (mutualinfo,corratio,normcorr,normmi,leastsq)
brain_seg.brain_mask_prob = 0.05; % Threshold from the prob maps to use to extract brain.
brain_seg.longSyN = 1; % Use long version of antsRegistrationSyN (way way longer, but better)
brain_seg.redebias = 0; % set to 1 to run debiasing on the extracted brain and re-run the first steps of brain extraction (might not be advised for heavy biased images)
brain_seg.rereg = 0; % set to 1 to re-run non-linear registration using the brain extraction obtained after segmentation & to re-do segmentation
brain_seg.use_brain_for_seg = 0; % set to 1 to perform the segmentation on the already extracted brain
brain_seg.use_new_segment = 1; % if set to 1 & if subject is human, spm_segment will be used to perform the segmentation (instead of spm_old_segment)


%% Parameters that are specific to realign_runs_FSL
realign_FSL.MC_method = 'mcflirt'; % Method used for motion correction: 'mcflirt' or 'ants'
realign_FSL.mcflirt_cost = 'normmi'; % cost function for mcflirt (mutualinfo, woods, corratio, normcorr, normmi, leastsquares)
realign_FSL.mcflirt_dof = 6; % degrees of freedom for mcflirt: 6 for rigid body (appropriate in most cases), 7 to compensate for global scale changes, 12 for all degrees of freedom
realign_FSL.IR_method = 'flirt'; % Method used for inter-run registration: 'flirt' or 'ants'
realign_FSL.IR_BET = 0; % When set to 1 averged runs will be BETed before inter-run registration with flirt
realign_FSL.IR_BETf = 0.5; % -f option of BET. (0->1); default=0.5; smaller values give larger brain outline estimates. 0.5 for human, 0.2 for macaques
realign_FSL.flirt_cost = 'normmi'; % cost function for inter-run registrationwith flirt (mutualinfo,corratio,normcorr,normmi,leastsq)
realign_FSL.flirt_dof = 6; % degrees of freedom for for inter-run registration with flirt: 6 for rigid body (appropriate in most cases), 7 to compensate for global scale changes, 12 for all degrees of freedom


%% Parameters that are specific to realign_runs_SPM
realign_flags.quality = 0.9; % Quality versus speed trade-off.  Highest quality (1)
realign_flags.sep = 2; % in voxels, the default separation to sample the images.
realign_flags.fwhm = 2; % in voxels, The FWHM applied to the images before estimating the realignment parameters.
realign_flags.interp = 2; % B-spline degree used for interpolation
realign_flags.wrap = [0 0 0]; % wrap around x, y or z. [0 0 0] for no wrap

%% Fieldmap correction
reslice_flags.unwarp = 1; % If set to 0, no fieldmap correction will be done

% SPM reslicing options (used by applyVDM in realign_runs_FSL, used by realign & unwarp or realign & reslice in realign_runs_SPM)
reslice_flags.mask = 1; % mask output images (true/false)
reslice_flags.mean = 0; % write mean image (true/false)
reslice_flags.interp = 4; % the B-spline interpolation method
reslice_flags.which = 2; % 2-reslice all the images
reslice_flags.prefix = 'resliced_'; % prefix for resliced images

%% Parameters of the algorithm that will classify each volume as steady or moving 
% Maximum translations allowed
mvt_params.max_trans_vox_x = 0.6; % in voxels, maximum translation allowed between 2 consecutive volumes (in real voxel size)
mvt_params.max_trans_vox_y = 0.6; % in voxels
mvt_params.max_trans_vox_z = 0.6; % in voxels
mvt_params.trans_fact = 2; % in numbers of max_trans_vox, maximum translation allowed between each volume and the mean of the steadiest period

% Maximum rotations allowed
mvt_params.max_rot_pitch = 0.25; % in degrees (0.3 seems ok), maximum rotation allowed between 2 consecutive volumes
mvt_params.max_rot_roll = 0.25; % in degrees
mvt_params.max_rot_yaw = 0.25; % in degrees
mvt_params.rot_fact = 2; % in numbers of max_rot, maximum rotation allowed between each volume and the mean of the steadiest period

% Steady period definition
mvt_params.min_dur_steady = 6; % in seconds, minimal duration of a steady period (shorter steady periods will be considered as movement)
mvt_params.min_dur_mvt = 1; % in seconds, minimal duration of a movement period (shorter movement periods will be considered as steady)
mvt_params.after_mvt_lag = 1; % in seconds, minimal duration to consider as a movement period after a detected movement (can be short for event-related protocol in which the stimulation is movement dependent. Should be around 10 sec for block designs)

% Verificaton scan
mvt_params.n_vols_verif = 3; % number of volumes per run to include in the verification scan

%% Smoothing parameters
smooth_params.fwhm = [2 2 2]; % in real voxels, Specify the FWHM of the Gaussian smoothing kernel in voxels. Three values should be entered, denoting the FWHM in the x, y and z directions
smooth_params.prefix = 'smoothed_'; %Specify the string to be prepended to the filenames of the smoothed image file(s). Default prefix is 's'.

%% Coregistration parameters
coreg_params.method = 'bbr'; % Method used to coregister the average functional volume to the anatomic volume:
% - 'bbr' will use FSL's epi_reg. This is the linear method that will give the best results, but it will only
%	 work with BOLD functional volumes (will most likely not work with MION volumes for instance) that are already very close
%	 to the anat (as is the case for human data). 
% - 'ants' will use the non-linear registration of ANTs. Suited for NHP with contrast agent (but also work with humans).

coreg_params.BETf = 0.5; % -f option of BET for functional average. (0->1); default=0.5; smaller values give larger brain outline estimates. 0.5 for human, 0.3 for macaques
coreg_params.flirt_cost = 'normmi'; % cost function for func-anat registration with flirt (mutualinfo,corratio,normcorr,normmi,leastsq)
coreg_params.use_close2func_anat = 0; % Set to 1 to use an anat scan already close to the func scans to help coregistration
coreg_params.c2f_BETf = 0.2; % -f option of BET for the anat that is close to the funcs. (0->1); default=0.5; smaller values give larger brain outline estimates. 0.5 for human, 0.3 for macaques



%% GLM parameters
GLM_params.mask_type = 'brain';  % Mask used when finding the t-threshold
							% '': no mask
							% 'brain': brain mask (created after brain segmentation)
							% 'grey': grey matter only (from the segmentation)
							
GLM_params.hrf = [5.4 5.2 10.8 7.35 0.35]; % The parameters of the hrf are specified by a row vector whose elements are:
									% 1. PEAK1: time to the peak of the first gamma density;
									% 2. FWHM1: approximate FWHM of the first gamma density;
									% 3. PEAK2: time to the peak of the second gamma density;
									% 4. FWHM2: approximate FWHM of the second gamma density;
									% 5. DIP: coefficient of the second gamma density.
									% The final hrf is: gamma1/max(gamma1)-DIP*gamma2/max(gamma2) scaled so that its total integral is 1.
									% If PEAK1=0 then there is no smoothing of that event type with the hrf.
									% If PEAK1>0 but FWHM1=0 then the design is lagged by PEAK1.
									% The default, chosen by Glover (1999) for an auditory stimulus, is: [5.4 5.2 10.8 7.35 0.35]
									% You can plot the HRF by running 'plot_hrf' after having ran this parameters file

GLM_params.pass_rejected_vols_as_regressors = 1; % set to 1 to create one regressor for each volume that was rejected during the realign step
GLM_params.all_trials_are_ok = 1; % set to 1 if each stimulus onset was only recorded when a trial was succesful (trial design).								
GLM_params.motion_reg_type = 0; % 0 = no motion regressors (MRs) ; 1 = 6 MRs (trans + rot) ; 2 = 12 MRs (6 + squared) ; 3 = 24 MRs (12 + derviate)
GLM_params.nPCs.white = 8; % Number of PCs from the PCA on white matter that will be passed as regressors
GLM_params.nPCs.csf = 9; % Number of PCs from the PCA on CSF that will be passed as regressors
GLM_params.GLMdenoise = 0; % use GLM_denoise to derive noise regressors instead of motion parameters
GLM_params.GLMdenoise_nPCs.white = 1; % number of PCs from the white matter PCA to use as extraregressors in GLMdenoise
GLM_params.GLMdenoise_nPCs.csf = 0; % number of PCs from the csf PCA to use as extraregressors in GLMdenoise
GLM_params.AR = 1; % the order of the autoregressive model. Order 1 seems to be adequate for 1.5T data, but higher order models may be needed for 3T data (will take much longer). 
GLM_params.n_temporal_trends = 3; % number of cubic spline temporal trends to be removed per 6 
								  %  minutes of scanner time (so it is backwards compatible). Temporal  
								  %  trends are modeled by cubic splines, so for a 6 minute run, N_TEMPORAL
								  %  <=3 will model a polynomial trend of degree N_TEMPORAL in frame times, 
								  %  and N_TEMPORAL>3 will add (N_TEMPORAL-3) equally spaced knots.
								  %  N_TEMPORAL=0 will model just the constant level and no temporal trends.
								  %  N_TEMPORAL=-1 will not remove anything, in which case the design matrix 
								  %  is completely determined by X_CACHE.X. 

p_val_peak = [0.05 0.01 0.001 0.0001 0.00001];
                                  
%% CONTRASTS
% The contrasts 'sound_vs_silence' as well as each condition vs silence are
% computed automatically
% Fill in the other contrasts that you wish to do here

contrasts.L = 2; % At which sound category level the analysis will be performed

% contrast names: better to name them cond1_vs_cond2 (no space allowed)
% contrast_names = {'voice_vs_nonvoice';...
% 				  'human_vs_nonvoice';...
%                   'macaque_vs_human';...
%                   'macaque_vs_all';...
%                   'human_vs_all';...
%                   'marmoset_vs_all'};
% contrast weights (the conditions order are displayed for each category level when running this file)
% contrast_weights = [ 1  1  1 -3  0;
% 					 1  0  0 -1  0;...
%                     -1  1  0  0  0;...
%                     -1  3 -1 -1  0;...
%                      3 -1 -1 -1  0;...
%                     -1 -1  3 -1  0];

% contrast names: better to name them cond1_vs_cond2 (no space allowed)
contrast_names = {'speech_vs_nonhuman';...
				  'nonspeech_vs_nonhuman';...
				  'speech_vs_nonspeech'};
% contrast weights (the conditions order are displayed for each category level when running this file)
contrast_weights = [ 6  6  0  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0;...
					 0  0  6  6 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0;...
					 1  1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0 0];


% % contrast names: better to name them cond1_vs_cond2 (no space allowed)
% contrast_names = {'male_vs_female'};
% % contrast weights (the conditions order are displayed for each category level when running this file)
% contrast_weights = [ 0 -1  1  0];


% % contrast names: better to name them cond1_vs_cond2 (no space allowed)
% contrast_names = {'artificial_vs_natural'};
% % contrast weights (the conditions order are displayed for each category level when running this file)
% contrast_weights = [0 1 0 0 -1 0 0 0 0 0 0 0 0];




%%% Restrict analysis to certain conditions
% contrasts.restrict(1).L = 1; % Level at which you want the restriction
% contrasts.restrict(1).c = [1 0 0 0 0]; % vector of zeros & ones, ones being the conditions at which the analysis will be restricted
% 
% contrasts.restrict(2).L = 3; % Add a new restriction by adding a new element (increase the number in brackets)
% contrasts.restrict(2).c = [0 0 0 0 0 0 1 1];



%%% Exclude certain conditions from the analysis
% contrasts.exclude(1).L = 1; % Level at which you want the exclusion
% contrasts.exclude(1).c = [0 0 1 0 0]; % vector of zeros & ones, ones being the excluded conditions
% 
% contrasts.exclude(2).L = 3; % Add a new exclusion by adding a new element (increase the number in brackets)
% contrasts.exclude(2).c = [0 0 0 0 0 0 1 1];






%% Parameters used to compute brain RDMs
paths.L = 2; % Category level at which you want the RDMs to be computed. You can comment or leave '' to compute RDMs at the stimulus level (or indicate the acutal level of the stim level).
% paths.CV_results_sub_folder = 'R2s64_group1'; % Results sub-folder where the contrast used for delimiting the clusters (ROIs) is. Use this for cross-validation, comment or leave '' otherwise.
RDM_params.contrast = 'human_vs_all'; % name of the contrast that will be used to mask the data
RDM_params.distance = 'euclidean'; % A string indicating the distance measure with which to calculate the RDMs (euclidean,correlation,spearman... [help pdist for other distances]).
RDM_params.peak_threshold = 2; % value at which this contrast t-map will thresholded
RDM_params.mask_with_grey = 0; % if set to 1, non-grey matter voxels will be removed from analysis
RDM_params.min_clust_size = 20; % only keep clusters from the contrast t-map that are larger than this value
RDM_params.n_voxels = 19; % size of the max cluster. Examples:
						  % 1: maximum only, 7: +faces, 19: +edges, 27: +corners
						  % 33: +faces2, 57: +edges2, 81: +corners2, 125: 5-cube
RDM_params.n_peaks = 0; % nb of maxima to keep in each cluster.
						  % if >= 1, same nb of max in each cluster
						  %	if < 1, the nb of peaks of each cluster will be:
						  % floor((sizeR / n_voxels) * n_peaks), where sizeR is the 
						  % size of the cluster being analyzed
						  % if == 0, all available maxima will be selected
RDM_params.colormap = 'parula'; % 'rsa' for rsa-toolbox colormap, or any Matlab colormap (parula,hot,jet...)





%% SCRIPT (do not modify)
addpath(genpath(fileparts(which('RunScript.m'))))
Pstack = dbstack;
set_parameters




%% Trash

% mvt_params.sel_from_fmap = 1; % Method to select the reference volume of each run:
%                               % '1' will select the volume, from the steadiest period, that has the minimum normmi cost when compared with the magnitude map (if available)
%                               % '0' will select the volume that is the closest to the average of the steadiest period during the run


