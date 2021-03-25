function out_file = xfm_temp2func(in_file,paths,coreg_params,out_file,interpf)

% Apply transforms to register to anat & template
setenv('ANTSPATH',paths.ANTS_path);
fprintf('Applying spatial transforms...\n\n')
if strcmp(coreg_params.method,'bbr')
    t_anat_to_average_warp = '';
elseif strcmp(coreg_params.method,'ants')
    t_anat_to_average_warp = [' -t ' paths.anat_to_average_warp];
else
    error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
end


if ~exist('out_file','var') || isempty(out_file)
    [~,in,ie] = fileparts(in_file);
    out_file = fullfile(paths.resliced,sprintf('In-%s-func_%s%s',paths.subject,in,ie));
end

if ~exist('interpf','var') || isempty(interpf)
	interpf = 'Linear';
end

system(sprintf('%s -i %s -r %s -o %s%s -t %s -t %s -t %s -n %s',...
    fullfile('$ANTSPATH','antsApplyTransforms'),...
    in_file,paths.average,out_file,...
    t_anat_to_average_warp,paths.anat_to_average_xfm,...
    paths.temp_to_anat_warp,paths.temp_to_anat_xfm,interpf));

fprintf('Output file:\n%s\n',out_file)            