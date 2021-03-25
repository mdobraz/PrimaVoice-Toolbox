% run parameters_... beforehand
%% Which parameters file should be used
this_params_file = Pstack(1).name;

%% Test each rank given by inter_session_cost

for max_rank = 1:n_total_runs
    fprintf('Computing costranks for %s & max_rank %i\n\n',this_params_file,max_rank)

    eval(this_params_file) % run parameters file

    run_GLM

    % clear & delete results_multi folder
    delete(paths.log_file_params)
    delete(sprintf('%s/*',paths.results_multi))
    rmdir(paths.results_multi)
end


