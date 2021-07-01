function out_file = xfm_func2anat(in_file,paths,coreg_params,out_file,interpf)


% Apply transforms to register to template
setenv('ANTSPATH',paths.ANTS_path);
fprintf('Applying spatial transforms...\n\n')
if strcmp(coreg_params.method,'bbr')
    t_average_to_anat_warp = '';
elseif strcmp(coreg_params.method,'ants')
    t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
else
    error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
end

if ~exist('out_file','var') || isempty(out_file)
    [~,in,ie] = fileparts(in_file);
    out_file = fullfile(paths.results_multi,sprintf('In-sub-%s-anat_%s%s',paths.subject,in,ie));
end

if ~exist('interpf','var') || isempty(interpf)
    interpf = 'Linear';
end

system(sprintf('%s -i %s -r %s -o %s%s -t %s -n Linear',...
    fullfile('$ANTSPATH','antsApplyTransforms'),...
    in_file,paths.anat.full,out_file,...
    t_average_to_anat_warp,paths.average_to_anat_xfm));


fprintf('Output file:\n%s\n',out_file)

% % Apply transforms to register to template
% in_file = sprintf('%s%s.nii.gz',out_base,ws{1}{st});
% out_file = fullfile(paths.results_multi,['In-' paths.template_name '_sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} ws{1}{st} '.nii.gz']);
% system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n Linear',...
%     fullfile('$ANTSPATH','antsApplyTransforms'),...
%     in_file,paths.anat.full,out_file,...
%     paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
%     t_average_to_anat_warp,paths.average_to_anat_xfm));


