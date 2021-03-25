% Launch this script from the head of frioul (not from a node), with a command
% like this one:
% frioul_batch -m "addpath('/hpc/banco/Primavoice_Scripts/Pipeline');RunScript"
% 'man frioul_batch' on frioul for more options

addpath(genpath(fullfile(fileparts(which('RunScript.m')),'functions')))


script_id = 'parameters_Maga19'; % this will appear in the subject of the emails. Usefull to distinguish several analyses

%% Email recipent address
% recipient = 'clementine.BODIN@univ-amu.fr'; % your email address
recipient = 'trapeau@gmail.com'; % your email address

%% Email transmission
[~,HN] = system('hostname');
HN = HN(1:end-1); % remove newline character at the end
sujet = sprintf('Script ''%s'' started on %s',script_id,HN);
ST = dbstack;
message = fileread(ST(1).file);
gmailme(sujet,message,recipient);

try
    tic;
    %% Scripts to run (comment / uncomment)
    % parameters_EloukBold
    % brain_segmentation
    parameters_Maga
    % realign_runs_SPM
    % realign_runs_FSL
    % smooth_runs
    % coregister_func2anat
    % run_GLM
    % launch_jackknife
    % run_GLM_by_trial
    % pref_maps

    % second_level_parameters_macaque
    brain_RDMs


    % compute_all_T2s

    %% Email transmission
    sujet = sprintf('Script ''%s'' finished on %s',script_id,HN);
    message = sprintf('after %1.0f seconds',toc);
    gmailme(sujet,message,recipient);
    exit
catch ME
    %% There was an error somewhere, email transmission
    sujet = sprintf('Error in Script ''%s'' launched on %s',script_id,HN);
    gmailme(sujet,getReport(ME,'extended','hyperlinks','off'),recipient);
    exit
end
