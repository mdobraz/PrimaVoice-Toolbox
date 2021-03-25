function select_jackrank(max_rank,this_params_file)

fprintf('Launching select_jackrank for %s & max_rank %i\n\n',this_params_file,max_rank)

eval(this_params_file) % run parameters file

run_GLM

% clear & delete results_multi folder
delete(paths.log_file_params)
delete(sprintf('%s/*',paths.results_multi))
rmdir(paths.results_multi)