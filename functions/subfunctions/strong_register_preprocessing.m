% function strong_register(in_file,ref_file,ref_brain_mask,method,paths,wmseg_file)
% brain_extraction script does not work when called from a function, so this is a script instead

[in_path,in_name,ext] = fileparts(in_file);
if strcmp(ext,'.gz'); [~,in_name] = fileparts(in_name); end

%% in_file preprocessing
switch method
	case 'autoBET' % this will only be done when the registration method chosen is 'bbr'
		% find center of gravity
    	cog = nan(1,3);
        for j = 1:3
            [~,cogchar] = system(sprintf('%sfslstats %s -C | cut -f%i -d'' ''',paths.FSL_prefix,in_file,j));
            cog(j) = round(str2double(cogchar));
        end
	
	    % Run BET
        BET_out_file = fullfile(in_path,[in_name '_BET.nii']);
	    if exist(BET_out_file,'file') == 2; delete(BET_out_file);end
	    system(sprintf('%sbet %s %s -c %i %i %i',paths.FSL_prefix,in_file,BET_out_file,cog(1),cog(2),cog(3)));
	    system(sprintf('gunzip %s',BET_out_file));
        
        view_slice_overlay(in_file,BET_out_file,0)
        figurewrite(fullfile(in_path,[in_name '_BET']),[],[],[],1) % Using GLMdenoise function
		
		params_file =  fullfile(in_path,[in_name '_BETparameters.json']);
	    jstruct.cog = cog;
	    spm_jsonwrite(params_file,jstruct)
	case {'manualBET','anat_manualBET'}
        
        if strcmp(method,'anat_manualBET') % crop volume if preprocessing an anat file
            force_type = [in_name ext];
            force_in_file = in_file;
            crop_volume
            fprintf('\n\nPress Enter when brain extraction is finished\n\n')
            pause
            in_file = fullfile(in_path,[in_name 'Cropped.nii']);
        end
        
		% N4BiasFieldCorrection
		[in_file,N4ed] = ANTS_N4(in_file,paths.ANTS_path);
		
		% Denoising
		[in_file,denoised] = spm_sanlm(in_file);
		
		% BET
		[in_path,in_name,ext] = fileparts(in_file);
		if strcmp(ext,'.gz'); [~,in_name,ext] = fileparts(in_name); end
		
		force_type = [in_name ext];
    	force_in_file = in_file;
	    force_out_path = in_path;
        BET_out_file = fullfile(in_path,[in_name '_BET.nii']);
	    brain_extraction
	    fprintf('\n\nPress Enter when brain extraction is finished\n\n')
        pause
	    clear force_out_path
	    
	otherwise
		error('''%s'' method does not exist. Available methods are:\n''autoBET'', ''manualBET'', ''anat_manualBET''',method)
end


%% TRASH

% % find center of gravity
%     	P = spm_vol(in_file);
%     	Y = spm_read_vols(P);
% 	    mask = true(size(Y));
% 	    measurements = regionprops(mask,Y,'WeightedCentroid');
% 	    cog = round(measurements.WeightedCentroid);

% system(sprintf('%sbet %s %s -c %i %i %i',paths.FSL_prefix,in_file,BET_out_file,cog(2),cog(1),cog(3))); % invert 1st & 2nd dims because "The first element of WeightedCentroid is the horizontal coordinate (or x-coordinate) of the weighted centroid. The second element is the vertical coordinate (or y-coordinate)."