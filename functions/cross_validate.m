group_folder = paths.results_multi(1:end-1);
a = dir([group_folder '*']);
n_groups = numel(a);
[~,group_name] = fileparts(group_folder);
CV_folder = fullfile(paths.results,['CV_' group_name 's']);


%% registration init
setenv('ANTSPATH',paths.ANTS_path);
if strcmp(coreg_params.method,'bbr')
    t_average_to_anat_warp = '';
elseif strcmp(coreg_params.method,'ants')
    t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
else
    error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
end


%% Contrasts
for con = 1:numel(contrasts.names)
    fprintf('%s:\n',contrasts.names{con})
    Yt = cell(n_groups,1);
    for g = 1:n_groups
        in_base = fullfile(sprintf('%s%i',group_folder,g),['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con}]);
        t_map_file = sprintf('%s_t.nii.gz',in_base);
        Pt = spm_vol(t_map_file);
        Yt{g} = spm_read_vols(Pt);
    end
    Xt = zeros(size(Yt{1}));

    % positive values
    idx = ones(size(Yt{1}));
    for g = 1:n_groups
        idx(Yt{g}<0) = 0; % reject negative values
    end
    idx = logical(idx);
    in_min = sprintf('Yt{%i}(idx),',1:n_groups);
    in_min(end) = '';
    eval(['Xt(idx) = min(' in_min ');']);

    % negative values
    idx = ones(size(Yt{1}));
    for g = 1:n_groups
        idx(Yt{g}>0) = 0; % reject positive values
    end
    idx = logical(idx);
    eval(['Xt(idx) = max(' in_min ');']);

    out_base = ['sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_t.nii'];

    Pt.fname = fullfile(CV_folder,out_base);
    % Pt = rmfield(Pt,'pinfo');
    delete([Pt.fname '*']);
    if ~exist(CV_folder,'dir');mkdir(CV_folder);end % create folder if non-existant
    spm_write_vol(Pt,Xt);

    system(['gzip ' Pt.fname]);
    % registration
    fprintf('Applying spatial transforms...\n\n')
    in_file = [Pt.fname '.gz'];
    out_file = fullfile(CV_folder,['In-' T1w_files.(paths.template_name).name '_sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_t.nii.gz']);
    system(sprintf('%s -i %s -r %s -o %s%s -t %s -n Linear',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        t_average_to_anat_warp,paths.average_to_anat_xfm));
 
    % maps to template
    out_file = fullfile(CV_folder,['In-' paths.template_name '_sub-' paths.subject '_res-' paths.results_name '_' contrasts.names{con} '_t.nii.gz']);
    system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n Linear',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
        t_average_to_anat_warp,paths.average_to_anat_xfm));
end



%% Preference maps
fprintf('\n\nPreference maps...\n\n')
if exist('Ypref','var')
    clear Ypref
    clear Yprob
end

for g = 1:n_groups
    pref_file = fullfile(sprintf('%s%i',group_folder,g),'Prefs.nii.gz');
    prob_file = fullfile(sprintf('%s%i',group_folder,g),'Probs.nii.gz');
    if ~exist(pref_file,'file')
        error('Cannot compute preference map intersection, no preference map was found.')
    end
    Ppref = spm_vol(pref_file);
    Pprob = spm_vol(prob_file);
    if ~exist('Ypref','var')
        Ypref = spm_read_vols(Ppref);
        Xpref = nan(size(Ypref));
        Yprob = spm_read_vols(Pprob);
    else
        Ypref(:,:,:,g) = spm_read_vols(Ppref);
        Yprob(:,:,:,g) = spm_read_vols(Pprob);
    end
end

for i = 1:size(Ypref,1)
    for j = 1:size(Ypref,2)
        for k = 1:size(Ypref,3)
            vect = squeeze(Ypref(i,j,k,:));
            if length(unique(vect)) > 1
                Xpref(i,j,k) = 0;
            else
                Xpref(i,j,k) = unique(vect);
            end
            Xprob = prod(Yprob,4);
        end
    end
end

Xprob(Xpref==0) = 0;

clear out_base
out_base{1} = 'Prefs_intersect.nii';
out_base{2} = 'Probs_intersect.nii';

Ppref.fname = fullfile(CV_folder,out_base{1});
delete([Ppref.fname '*']);
spm_write_vol(Ppref,Xpref);


Pprob.fname = fullfile(CV_folder,out_base{2});
delete([Pprob.fname '*']);
spm_write_vol(Pprob,Xprob);


%% Registration
fprintf('Applying spatial transforms...\n\n')
interp = {'NearestNeighbor';'Linear'};
for i = 1:2
    system(['gzip ' fullfile(CV_folder,out_base{i})]);
    % Apply transforms to register to anat & template
    in_file = fullfile(CV_folder,[out_base{i} '.gz']);
    out_file = fullfile(CV_folder,['In-' T1w_files.(paths.template_name).name '_' out_base{i} '.gz']);
    system(sprintf('%s -i %s -r %s -o %s%s -t %s -n %s',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        t_average_to_anat_warp,paths.average_to_anat_xfm,...
        interp{i}));


    out_file = fullfile(CV_folder,['In-' paths.template_name '_' out_base{i} '.gz']);
    system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n %s',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
        t_average_to_anat_warp,paths.average_to_anat_xfm,...
        interp{i}));
end