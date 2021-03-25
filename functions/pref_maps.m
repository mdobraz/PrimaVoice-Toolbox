%%%%% Compute preference map %%%%%
fprintf('\n### Preference maps ###\n')
nconds = contrasts.nconds;

%% Get dims
P = spm_vol(paths.average); % just to get vol dim & a P
dims = P.dim;
P = rmfield(P,'pinfo');

%% Compute prefs
conds = 2:nconds;
[Prefs,Probs] = compute_prefs(conds,contrasts,dims,paths);

% write images
clear out_base
out_base{1} = 'Prefs.nii';
out_base{2} = 'Probs.nii';

P.fname = fullfile(paths.results_multi,out_base{1});
P.descrip = 'General preferences from contrasts vs_silence';
delete([P.fname '*']);
spm_write_vol(P,Prefs);


P.fname = fullfile(paths.results_multi,out_base{2});
P.descrip = 'Probabilities of significant preference from contrasts vs_silence';
delete([P.fname '*']);
spm_write_vol(P,Probs);


%% Registration
setenv('ANTSPATH',paths.ANTS_path);
fprintf('Applying spatial transforms...\n\n')
if strcmp(coreg_params.method,'bbr')
    t_average_to_anat_warp = '';
elseif strcmp(coreg_params.method,'ants')
    t_average_to_anat_warp = [' -t ' paths.average_to_anat_warp];
else
    error('''coreg_params.method'' must be ''bbr'' or ''ants''. Modify your parameters file')
end
interp = {'NearestNeighbor';'Linear'};
for i = 1:2
    system(['gzip ' fullfile(paths.results_multi,out_base{i})]);
    % Apply transforms to register to anat & template
    in_file = fullfile(paths.results_multi,[out_base{i} '.gz']);
    out_file = fullfile(paths.results_multi,['In-' T1w_files.(paths.template_name).name '_' out_base{i} '.gz']);
    system(sprintf('%s -i %s -r %s -o %s%s -t %s -n %s',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        t_average_to_anat_warp,paths.average_to_anat_xfm,...
        interp{i}));


    out_file = fullfile(paths.results_multi,['In-' paths.template_name '_' out_base{i} '.gz']);
    system(sprintf('%s -i %s -r %s -o %s -t %s -t %s%s -t %s -n %s',...
        fullfile('$ANTSPATH','antsApplyTransforms'),...
        in_file,paths.anat.full,out_file,...
        paths.anat_to_temp_warp,paths.anat_to_temp_xfm,...
        t_average_to_anat_warp,paths.average_to_anat_xfm,...
        interp{i}));
end

return










%% Compute preference maps corresponding to each contrast
pref_cons = nconds+1:length(con_weights); % will determine preferences based on the last contrasts
for icon = 1:length(pref_cons)
    con = pref_cons(icon);
    
    fprintf('\t%s...',con_names{con})
    
    ind = con_weights{con}<0;
    ind = ind(2:end);
    if sum(ind) == 1
        weak_con = silence_Y(:,:,:,ind);
    elseif sum(ind) == 2
        weak_con = max(silence_Y(:,:,:,ind),[],4);
    else
        error('Cannot compute preferences for this contrast')
    end
    weak_con_vect = weak_con(~isnan(weak_con));
    
    ind = con_weights{con}>0;
    ind = ind(2:end);
    if sum(ind) == 1
        strong_con = silence_Y(:,:,:,ind);
    elseif sum(ind) == 2
        strong_con = min(silence_Y(:,:,:,ind),[],4);
    else
        error('Cannot compute preferences for this contrast')
    end
    strong_con_vect = strong_con(~isnan(strong_con));


    pref = zeros(size(strong_con,1),size(strong_con,2),size(strong_con,3));
    for i = 1:size(strong_con,1)
        for j = 1:size(strong_con,2)
            for k = 1:size(strong_con,3)
                if ~isnan(strong_con(i,j,k))
                    diff_SW = strong_con(i,j,k) - weak_con(i,j,k);
                    [~,I] = min(abs(xpd - diff_SW));
                    pref(i,j,k) = csypd(I);
                end
            end
        end
    end

    % write image
    pref_file = sprintf('%s/prefs_%s.nii',paths.results,con_names{con});

    P.fname = pref_file;
    P.descrip = con_names{con};
    spm_write_vol(P,pref);
    fprintf(' done\n')
end


for i = 1:numel(out_files)
    [~,out_name] = fileparts(out_files{i});
    out_file = fullfile(paths.results_multi,[out_name '.nii']);
    if exist(out_file,'file') == 2; delete(out_file);end
    if coreg_params.use_daily_anat && ismember(coreg_params.daily_anat_method,{'fnirt','fnirt_brains'})
        system(sprintf('%sapplywarp --ref=%s --in=%s --out=%s --warp=%s',paths.FSL_prefix,paths.anat_file,out_files{i},out_file,paths.average_to_anat_warp));
    else
        system(sprintf('%sflirt -in %s -ref %s -out %s -applyxfm -init %s',paths.FSL_prefix,out_files{i},paths.anat_file,out_file,paths.average_to_anat_xfm));
    end
end

    
    
    
    


