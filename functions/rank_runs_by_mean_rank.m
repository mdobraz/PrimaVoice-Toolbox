%% MCxR2s
max_rank = input('R2s_Run_selection rank: ','s');
filename = fullfile(paths.subject_analysis,['R2s_Run_selection_' max_rank '.mat']);
R2s = load(filename);

max_rank = input('MotionCorr_Run_selection rank: ','s');
filename = fullfile(paths.analysis,['MotionCorr_Run_selection_' max_rank '.mat']);
MC = load(filename);

if length(R2s.classement) == length(MC.classement)
	[selected_sessions,selected_runs,classement,max_rank,selected] = plot_rank_corr(R2s,MC,'Relaxation rate ranking','Motion RMS ranking');

	% Save
	fig_prefix = fullfile(paths.analysis,['MCxR2s_Run_selection_' num2str(max_rank)]);
	MCxR2s_filename = [fig_prefix '.mat'];
	save(MCxR2s_filename,'selected_sessions','selected_runs','classement','max_rank')
	fprintf('\nSelection saved in:\n%s\n\n',MCxR2s_filename)
	input('### Press enter to close figure and finish script')
	figurewrite(fig_prefix,[],0,paths.analysis); % the 0 is to force eps figure
else
	warning('R2s and Motion run selection were not done on the same number of runs\n%i runs for R2s\n%i runs for Motion',length(R2s.classement),length(MC.classement))
end



%% Comparison with Jackknife
max_rank = input('Jackknife_Run_selection rank: ','s');
filename = fullfile(paths.analysis,['Jackknife_Run_selection_' max_rank '.mat']);
if ~exist(filename,'file')
	error('Jackknife run selection file does not exist.')
end
jack = load(filename);

if exist('MCxR2s_filename','var') && exist(MCxR2s_filename,'file')
	MCxR2s = load(MCxR2s_filename);
	plot_rank_corr(jack,MCxR2s,'Jackknife ranking','MC x R2s',selected);
end

if length(jack.classement) == length(MC.classement)
	plot_rank_corr(jack,MC,'Jackknife ranking','Motion RMS ranking');
else
	warning('Jackknife and Motion run selection were not done on the same number of runs\n%i runs for Jackknife\n%i runs for Motion',length(jack.classement),length(MC.classement))
end

if length(jack.classement) == length(R2s.classement)
	plot_rank_corr(jack,R2s,'Jackknife ranking','Relaxation rate ranking');
else
	warning('Jackknife and R2s run selection were not done on the same number of runs\n%i runs for Jackknife\n%i runs for R2s',length(jack.classement),length(R2s.classement))
end



