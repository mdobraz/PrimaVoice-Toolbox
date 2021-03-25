clear all
sourcedata = '/hpc/banco/Primavoice_Data_and_Analysis/sourcedata';

newset_char = tdfread(fullfile(sourcedata,'TESTING-NEW.csv'));
fn = fieldnames(newset_char);

for i = 1:length(fn)
	for j = 1:size(newset_char.(fn{i}),1)
		newset.(fn{i}){j,1} = strtrim(newset_char.(fn{i})(j,:));
	end
end


a = dir(fullfile(sourcedata,'sub-*'));

for sub = 1:length(a)
	subdir = fullfile(sourcedata,a(sub).name)
	b = dir(fullfile(subdir,'ses-*'));
	for s = 1:length(b)
		sesdir = fullfile(subdir,b(s).name,'func')
		c = dir(fullfile(sesdir,'*ScanLog.mat'));
		for f = 1:length(c)
			oldfile = fullfile(sesdir,c(f).name);
			newfile = fullfile(sesdir,['orig_' c(f).name]);
			copyfile(oldfile,newfile);
			load(oldfile)

			for i = 1:length(newset.L5)
				ScanLog.fMRISTAT.L5.names{i} = newset.L5{i};
				ScanLog.SPM.L5.names{i} = newset.L5{i};
			end
			uL1 = unique(newset.L1);
			for i = 1:length(uL1)
				ScanLog.fMRISTAT.L1.names{i} = uL1{i};
				ScanLog.SPM.L1.names{i} = uL1{i};
			end
			uL2 = unique(newset.L2);
			uL2{end+1} = 'silence';
			ScanLog.fMRISTAT.L2.names = uL2;
			ScanLog.SPM.L2.names = uL2;
			ScanLog.SPM.L2.onsets = cell(length(uL2),1);
			ScanLog.SPM.L2.durations = cell(length(uL2),1);
			for i = 1:size(ScanLog.fMRISTAT.L5.events,1)
				sn = ScanLog.fMRISTAT.L5.events(i,1);
				if sn > length(newset.L2)
					L2name = 'silence';
				else
					L2name = newset.L2{sn};
				end
				indexC = strfind(uL2,L2name);
				index = find(not(cellfun('isempty',indexC)));
				ScanLog.fMRISTAT.L2.events(i,1) = index;
				ScanLog.SPM.L2.sequence(i) = index;
			end

			save(oldfile,'ScanLog')
		end
	end
end