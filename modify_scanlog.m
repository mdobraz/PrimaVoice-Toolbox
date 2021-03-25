paths.subject = 'EloukMION';
paths.task = 'BlocLocalizer';
% sessions = [16 18:26];
sessions = 2:6;

data_folder = ['/hpc/banco/Primavoice_Data_and_Analysis/sub-' paths.subject];
paths.tsv_sequences = ['/hpc/banco/Primavoice_Data_and_Analysis/sourcedata/sub-' paths.subject];


for s = 1:length(sessions)
	session = sessions(s);
	ses_folder = fullfile(data_folder,sprintf('ses-%02.0f',session),'func');
	a = dir(fullfile(ses_folder,'*.json'));
	for i = 1:numel(a)
		[~,name] = fileparts(a(i).name);
		scan_info = spm_jsonread(fullfile(ses_folder,[name '.json']));
		run = str2num(name(strfind(name,'_run-')+5:strfind(name,'_run-')+6));
		create_conditions_matrix_no_mvt(paths,session,run,scan_info)
	end
end