%% Setting paths & -p option for bash scripts
debias_path = which('T1xT2BiasFieldCorrection.sh');
BET_path = which('T1xT2BET.sh');
IterREGBET_path = which('IterREGBET.sh');
c3d_path = which('c3d_affine_tool.b');
if isfield(paths,'FSL_prefix') && ~isempty(paths.FSL_prefix)
    popt = ['-p ' paths.FSL_prefix];
else
    popt ='';
end


%% Close 2 func anat
if coreg_params.use_close2func_anat
	if strcmp(coreg_params.method,'bbr')
		error('Usage of intermediary anat is not compatible with the ''bbr'' option. Change your parameters file.')
	end
	if ~isfield(paths,'close2func_anat_session') || isempty(paths.close2func_anat_session)
		load(paths.reference_scan_infos)
		c2f_session = Master.session;
	else
		c2f_session = paths.close2func_anat_session;
	end
	c2f_json = get_anat(BIDS,paths,paths.close2func_anat_suffix,c2f_session);
	c2f_file = spm_jsonread(c2f_json);

	c2f_path = fullfile(paths.segmentation,c2f_file.name); % folder where the preprocessing of the anat will happen
	if ~exist(c2f_path,'dir');mkdir(c2f_path);end % create folder if non-existant

	%% Denoising
	T1_denoised = spm_sanlm(c2f_file.raw,c2f_path);

	%% First large Debias (brain_seg.BETf * 0.6)
	system(sprintf('bash %s -t1 %s -t2 %s -bet 2 -f %s %s',debias_path,T1_denoised,T1_denoised,num2str(coreg_params.c2f_BETf * 0.6),popt));
	T1_debiased = [fileparts2(T1_denoised) '_debiased.nii.gz'];
	T1_BET_mask = [fileparts2(T1_debiased) '_BET_mask.nii.gz'];
	view_slice_overlay(T1_debiased,T1_BET_mask,0,[],[],[],lines(2))
	figurewrite(fileparts2(T1_BET_mask),[],[],[],1) % Using GLMdenoise function

	%% T1xT2BET
	system(sprintf('bash %s -t1 %s -t2 %s -n 3 -f %s %s',BET_path,T1_debiased,T1_debiased,num2str(coreg_params.c2f_BETf),popt));
	T1_BET = [fileparts2(T1_debiased) '_BET.nii.gz'];
	rel_cog = vol_rel_cog(T1_BET,paths);
	view_slice_overlay(T1_debiased,T1_BET,0,[],[],[],[],rel_cog)
	figurewrite(fileparts2(T1_BET),[],[],[],1) % Using GLMdenoise function

	%% IterREGBET
	system(sprintf('bash %s -inw %s -inb %s -refb %s -cost %s -dof 6 %s',IterREGBET_path,T1_debiased,T1_BET,paths.anat.brain,coreg_params.flirt_cost,popt));
    T1_brain = [fileparts2(T1_debiased) '_IRbrain.nii.gz'];
    rel_cog = vol_rel_cog(T1_brain,paths);
    view_slice_overlay(T1_debiased,T1_brain,0,[],[],[],[],rel_cog)
    figurewrite(fileparts2(T1_brain),[],[],[],1) % Using GLMdenoise function

    IR_ref = T1_brain;
else
	IR_ref = paths.anat.brain;
end



%% Apply bias field correction to the grand Average
system(sprintf('bash %s -t1 %s -t2 %s -bet 2 -f %s %s',debias_path,paths.average,paths.average,num2str(coreg_params.BETf * 0.6),popt));
avg_debiased = [fileparts2(paths.average) '_debiased.nii.gz'];
avg_BET_mask = [fileparts2(avg_debiased) '_BET_mask.nii.gz'];
view_slice_overlay(avg_debiased,avg_BET_mask,0,[],[],[],lines(2))
figurewrite(fileparts2(avg_BET_mask),[],[],[],1) % Using GLMdenoise function


%% BETing debiased averages, T1xT2BET
system(sprintf('bash %s -t1 %s -t2 %s -n 3 -f %f %s',BET_path,avg_debiased,avg_debiased,coreg_params.BETf,popt));
avg_BET = [fileparts2(avg_debiased) '_BET.nii.gz'];
rel_cog = vol_rel_cog(avg_BET,paths);
view_slice_overlay(avg_debiased,avg_BET,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(avg_BET),[],[],[],1) % Using GLMdenoise function


%% Apply anat brain mask to debiased averages (IterREGBET)
system(sprintf('bash %s -inw %s -inb %s -refb %s -dof 6 %s',IterREGBET_path,avg_debiased,avg_BET,IR_ref,popt));
avg_brain = [fileparts2(avg_debiased) '_IRbrain.nii.gz'];
rel_cog = vol_rel_cog(avg_brain,paths);
view_slice_overlay(avg_debiased,avg_brain,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(avg_brain),[],[],[],1) % Using GLMdenoise function




setenv('ANTSPATH',paths.ANTS_path);

if strcmp(coreg_params.method,'bbr')
	%% BBR
	system(sprintf('%sepi_reg --wmseg=%s --epi=%s --t1=%s --t1brain=%s --out=%s',...
		paths.FSL_prefix,paths.anat.tissues.white,...
		avg_brain,paths.anat.full,paths.anat.brain,paths.average_to_anat));
	[~,xfm_i2r_file] = fileparts(fileparts2(paths.average_to_anat));
	xfm_i2r_file = fullfile(paths.resliced,[xfm_i2r_file '.mat']);
	% movefile(xfm_i2r_file,paths.average_to_anat_xfm)


	%% Convert bbr affine to ants ITK format & invert it
	system(sprintf('%s -ref %s -src %s %s -fsl2ras -oitk %s',c3d_path,...
		paths.anat.brain,avg_brain,xfm_i2r_file,paths.average_to_anat_xfm));

	system(sprintf('%s -itk %s -inv -oitk %s',c3d_path,paths.average_to_anat_xfm,paths.anat_to_average_xfm));
	t_anat_to_average_warp = '';
	t_average_to_anat_warp = '';



elseif strcmp(coreg_params.method,'ants')


	%% Registration names (from IterREGBET)
	[~,anat_brain_name,ext] = fileparts(paths.anat.brain);
	if strcmp(ext,'.gz'); [~,anat_brain_name] = fileparts(anat_brain_name); end


	[~,avg_name] = fileparts(paths.average);
	debiased_avg_name = [avg_name '_debiased'];
	BETed_avg_name = [debiased_avg_name '_BET'];

	[~,avg_brain_name] = fileparts(avg_brain);
	[~,avg_brain_name] = fileparts(avg_brain_name);


	avg2anat_xfm = fullfile(paths.resliced,[BETed_avg_name '_FLIRT-to_' anat_brain_name '.xfm']);
	% anat2avg_xfm = fullfile(paths.resliced,[BETed_avg_name '_FLIRT-to_' anat_brain_name '_inverse.xfm']);

	if coreg_params.use_close2func_anat
		[~,BETed_c2f_name] = fileparts(T1_BET);
		[~,BETed_c2f_name] = fileparts(BETed_c2f_name);
		c2f2anat_xfm = fullfile(c2f_path,[BETed_c2f_name '_FLIRT-to_' anat_brain_name '.xfm']);
		[~,brain_c2f_name] = fileparts(T1_brain);
		[~,brain_c2f_name] = fileparts(brain_c2f_name);
		avg2c2f_xfm = fullfile(paths.resliced,[BETed_avg_name '_FLIRT-to_' brain_c2f_name '.xfm']);
		system(sprintf('%sconvert_xfm -omat %s -concat %s %s',paths.FSL_prefix,avg2anat_xfm,c2f2anat_xfm,avg2c2f_xfm));
	end


	avg_brain_in_anat = fullfile(paths.resliced,[BETed_avg_name '_FLIRT-to_' anat_brain_name '.nii.gz']);
	avg_in_anat = fullfile(paths.resliced,[debiased_avg_name '_FLIRT-to_' anat_brain_name '.nii.gz']);

	%% Apply transform to whole head average (func2anat)
	system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s',paths.FSL_prefix,avg_debiased,paths.anat.brain,avg_in_anat,avg2anat_xfm));
	if coreg_params.use_close2func_anat
		system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s',paths.FSL_prefix,avg_brain,paths.anat.brain,avg_brain_in_anat,avg2anat_xfm));
	end

	%% Run antsRegistrionSyNQuick.sh to find non-linear warp between registered avg & anat
	ANTs_out_base = fullfile(paths.resliced,[avg_brain_name '_SyN-to_' anat_brain_name '_']);
	system(sprintf('%s -d 3 -f %s -f %s -m %s -m %s -o %s -j 1',fullfile(paths.ANTS_scripts,'antsRegistrationSyNQuick.sh'),paths.anat.brain,paths.anat.full,avg_brain_in_anat,avg_in_anat,ANTs_out_base));


	%% Convert IterREGBET affine to ants ITK format
	avg2anat_mat = fullfile(paths.resliced,[BETed_avg_name '_AFFINE-to_' anat_brain_name '.mat']);
	system(sprintf('%s -ref %s -src %s %s -fsl2ras -oitk %s',c3d_path,paths.anat.brain,avg_BET,avg2anat_xfm,avg2anat_mat));

	%% Multiply FSL & ANTS affine matrices & invert it
	ants_affine = [ANTs_out_base '0GenericAffine.mat'];
	system(sprintf('%s -itk %s -itk %s -mult -oitk %s',c3d_path,avg2anat_mat,ants_affine,paths.average_to_anat_xfm));
	system(sprintf('%s -itk %s -inv -oitk %s',c3d_path,paths.average_to_anat_xfm,paths.anat_to_average_xfm));


	%% Move files
	movefile([ANTs_out_base '1Warp.nii.gz'],paths.average_to_anat_warp)
	movefile([ANTs_out_base '1InverseWarp.nii.gz'],paths.anat_to_average_warp)
	movefile([ANTs_out_base 'Warped.nii.gz'],paths.average_to_anat)

	t_anat_to_average_warp = [' -t ' paths.anat_to_average_warp];
	t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
else
	error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
end




%% Apply transforms to PCA masks
ref_file = avg_debiased;
pca_tissues = {'white';'csf'};
for i = 1:numel(pca_tissues)
    in_file = paths.anat.pca_mask.(pca_tissues{i});
    out_file = paths.func.pca_mask.(pca_tissues{i});
    system(sprintf('%s -i %s -r %s -o %s%s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),...
    	in_file,ref_file,out_file,t_anat_to_average_warp,paths.anat_to_average_xfm));
end

%% Apply transforms to brain mask
in_file = paths.anat.brain_mask;
out_file = paths.func.brain_mask;
system(sprintf('%s -i %s -r %s -o %s%s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),...
	in_file,ref_file,out_file,t_anat_to_average_warp,paths.anat_to_average_xfm));



%% Apply transforms to each tissue
for i = 1:numel(tissues)
    in_file = paths.anat.tissues.(tissues{i});
    out_file = paths.func.tissues.(tissues{i});
    system(sprintf('%s -i %s -r %s -o %s%s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),...
    	in_file,ref_file,out_file,t_anat_to_average_warp,paths.anat_to_average_xfm));
    
    probs = paths.anat.tissue_probs;
    for p = 1:numel(probs)
        prob = probs{p};
        in_file = paths.anat.tissues.(sprintf('p%s',probs{p})).(tissues{i});
        out_file = paths.func.tissues.(sprintf('p%s',probs{p})).(tissues{i});
        system(sprintf('%s -i %s -r %s -o %s%s -t %s -n NearestNeighbor',fullfile('$ANTSPATH','antsApplyTransforms'),...
    	in_file,ref_file,out_file,t_anat_to_average_warp,paths.anat_to_average_xfm));
    end
end


%% Apply transforms to template
ref_file = paths.template_brain;
in_file = avg_brain;
out_file = paths.average_to_temp;
system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n Linear',fullfile('$ANTSPATH','antsApplyTransforms'),...
	in_file,ref_file,out_file,paths.anat_to_temp_warp,paths.anat_to_temp_xfm,t_average_to_anat_warp,paths.average_to_anat_xfm));

% system(sprintf('%s -i %s -r %s -o %s%s -t %s -t %s -t %s -n Linear',fullfile('$ANTSPATH','antsApplyTransforms'),in_file,ref_file,out_file,t_average_to_anat_warp,paths.average_to_anat_xfm,paths.anat_to_temp_warp,paths.anat_to_temp_xfm));



%% Create full filled volume (mask used when no mask is specified)
if exist(paths.no_mask,'file') == 2; delete(paths.no_mask);end
system(sprintf('%sfslmaths %s -thr 0 %s -odt input',paths.FSL_prefix,paths.average,paths.no_mask));
system(sprintf('gunzip %s',paths.no_mask));


%% Check figure
view_slice_overlay(paths.template_head,paths.average_to_temp,0)
figurewrite(fileparts2(paths.average_to_temp),[],[],[],1) % Using GLMdenoise function

view_slice_overlay(avg_debiased,paths.func.pca_mask.white,0,[],[],[],lines(2),rel_cog)
figurewrite(fileparts2(paths.func.pca_mask.white),[],[],[],1) % Using GLMdenoise function

rel_cog = vol_rel_cog(paths.anat.brain,paths);
rel_cog(1) = rel_cog(1)*1.5;
view_slice_overlay(paths.anat.full,paths.average_to_anat,0,[],[],[],[],rel_cog)
figurewrite(fileparts2(paths.average_to_anat),[],[],[],1) % Using GLMdenoise function



%% Move figures
if ~exist(paths.resliced_fig,'dir');mkdir(paths.resliced_fig);end % create folder if non-existant
movefile(fullfile(paths.resliced,'*.png'),paths.resliced_fig);




%% Trash


