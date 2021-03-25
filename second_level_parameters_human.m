clearvars('-except','script_id','recipient','HN','run_index','max_rank','this_params_file'); % clear all but 'recipient' and 'HN' which are used in RunScript & 'run_index' used for jackknife procedure

%% PARAMETERS
% All the parameters of the analysis are set here
% Run this script before starting any analysis or using any function. Re-run it if you made changes.

%% Paths
paths.dataset = '/hpc/banco/Primavoice_Data_and_Analysis'; % Directory where the data are (general BIDS folder containing all subjects)
paths.subjects = {'02';'04';'05';'06';'07'}; % cell array of the subjects you want to analyze (eg. if the folder containing the subject data is called 'sub-moumou', omit the prefix 'sub-' and leave only 'moumou')
paths.task = 'sparse'; % the task you want to analyze (eg. if the datafiles include 'task-sparse', omit the prefix 'task-' and leave only 'sparse')
paths.analysis_name = 'human'; % name used to differentiate several second level analyses

%% Segmentation template (for non-human primates). Comment for humans.
paths.templates = '/hpc/banco/Primavoice_Data_and_Analysis/templates'; % folder where templates are stored
paths.template_name = 'MNI152'; % Template that will be used during the segmentation process. Name should match a 'case' in the switch found in functions/templates.m
                             % Edit file functions/templates.m to add a new template
%% FSL & ANTS paths
paths.FSL_prefix = 'fsl5.0-'; % prefix to call FSL tools by the system (eg. on frioul, to call 'flirt' you have to type 'fsl5.0-flirt', so the prefix is 'fsl5.0-'). Leave '' if no prefix is needed
paths.ANTS_path = '/hpc/soft/ANTS/antsbin/bin'; % path to ANTS binaries
paths.ANTS_scripts = '/hpc/soft/ANTS/ANTs/Scripts'; % path to ANTS scripts
   


%% Names of the contrasts at which the second level analysis will be performed
% comment to perform the analysis at the condition level or to compute brain RDMs

% contrast_names = {'sound_vs_silence';...
% 				  'voice_vs_nonvoice';...
% 				  'human_vs_nonvoice';...
%                   'macaque_vs_human';...
%                   'macaque_vs_all';...
%                   'human_vs_all';...
%                   'marmoset_vs_all'};


% contrast_names = {'human_vs_silence';...
% 				  'macaque_vs_silence';...
% 				  'marmoset_vs_silence';...
%                   'nonvocal_vs_silence'};

% contrast_names = {'speech_vs_nonhuman';...
% 				  'nonspeech_vs_nonhuman';...
% 				  'speech_vs_nonspeech'};

contrasts.L = 2; % Level at which the analysis will be performed if condition level is chosen (contrast_names commented out)


%% Parameters used to compute brain RDMs
paths.L = 2; % Category level at which you want the RDMs to be computed. You can comment or leave '' to compute RDMs at the stimulus level (or indicate the acutal level of the stim level).
% paths.CV_results_sub_folder = 'human_group1'; % Results sub-folder where the contrast used for delimiting the clusters (ROIs) is. Use this for cross-validation, comment or leave '' otherwise.
RDM_params.contrast = 'human_vs_all'; % name of the contrast that will be used to mask the data
RDM_params.distance = 'euclidean'; % A string indicating the distance measure with which to calculate the RDMs (euclidean,correlation,spearman... [help pdist for other distances]).
RDM_params.peak_threshold = 3.09; % value at which this contrast t-map will thresholded
RDM_params.mask_with_grey = 0; % if set to 1, non-grey matter voxels will be removed from analysis
RDM_params.min_clust_size = 300; % only keep clusters from the contrast t-map that are larger than this value
RDM_params.n_voxels = 297; % size of the max cluster. Examples:
						  % 1: maximum only, 7: +faces, 19: +edges, 27: +corners
						  % 33: +faces2, 57: +edges2, 81: +corners2, 125: 5-cube
RDM_params.n_peaks = 0; % nb of maxima to keep in each cluster.
						  % if >= 1, same nb of max in each cluster
						  %	if < 1, the nb of peaks of each cluster will be:
						  % floor((sizeR / n_voxels) * n_peaks), where sizeR is the 
						  % size of the cluster being analyzed
						  % if == 0, all available maxima will be selected
RDM_params.colormap = 'parula'; % 'rsa' for rsa-toolbox colormap, or any Matlab colormap (parula,hot,jet...)
RDM_params.fixed_random = 'fixed'; % char, use 'fixed' effects or 'random' effects group analysis





% plot_brain_RDMs(RDM_params.contrast,RDM_params.peak_threshold,RDM_params.mask_with_grey,RDM_params.min_clust_size,RDM_params.n_peaks,RDM_params.n_voxels,paths,ScanLog,RDM_params.fixed_random)







%% SCRIPT (do not modify)
addpath(genpath(fileparts(which('RunScript.m'))))
Pstack = dbstack;
set_second_level_parameters
