% run parameters_... beforehand

%% Which parameters file should be used
this_params_file = Pstack(1).name;

%% Running jackknife
for run_index = 1:n_total_runs
    fprintf('Launching jackknife for %s & run index %i\n\n',this_params_file,run_index)
    eval(this_params_file) % run parameters file

    run_GLM

    % clear & delete results_multi folder
    delete(paths.log_file_params)
    delete(sprintf('%s/*',paths.results_multi))
    rmdir(paths.results_multi)
end


%% Test each rank given by jackknife
jackknife_results
for max_rank = 1:numjack
    select_jackrank(max_rank,this_params_file)
end







%% TRASH

% function launch_jackknife(run_index)
% 
% % If you want to launch it from frioul (and have multiple jobs at the same time)
% % frioul_batch -c 1 -M "[r_[1:49]]" -m "addpath('/hpc/banco/Primavoice_Scripts/Pipeline');launch_jackknife"
% % "[r_[1:49]]" upper bound should be n_runs + 1
% 
% 
% fprintf('Launching jackknife for %s & run index %i\n\n',ST(1).name,run_index)
% 
% eval(ST(1).name)
% parameters_Maga
% 
% run_GLM
% 
% % clear & delete results_multi folder
% delete(paths.log_file_params)
% delete(sprintf('%s/*',paths.results_multi))
% rmdir(paths.results_multi)
% 
% 
% exit % if launched from frioul (comment otherwise)




% function select_jackrank(max_rank)
% 
% % If you want to launch it from frioul (and have multiple jobs at the same time)
% % frioul_batch -c 1 -M "[r_[1:49]]" -m "addpath('/hpc/banco/Primavoice_Scripts/Pipeline');select_jackrank"
% % "[r_[1:49]]" upper bound should be n_runs + 1
% 
% fprintf('Launching select_jackrank for max_rank %i\n\n',max_rank)
% 
% parameters_Maga
% 
% run_GLM
% 
% % clear & delete results_multi folder
% delete(paths.log_file_params)
% delete(sprintf('%s/*',paths.results_multi))
% rmdir(paths.results_multi)
% 
% 
% % exit % if launched from frioul (comment otherwise)

