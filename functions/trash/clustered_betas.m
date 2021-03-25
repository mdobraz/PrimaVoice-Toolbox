%%%%% Compute preference map %%%%%
fprintf('\n### Clustered betas ###\n')
nconds = contrasts.nconds;

%% Get dims
P = spm_vol(paths.average); % just to get vol dim & a P
dims = P.dim;

%% Pref contrasts
pref_cons = nconds+nconds+1:length(contrasts.names); % will determine preferences based on the last contrasts
pref_cons  = [1 pref_cons ]; % add sound_vs_silence contrast

% pref_cons = [1 12 14];

%% Compute prefs vs silence
conds = 2:nconds;
% conds = [2 3 5];
% compute_clustered_betas(p_val_peak(1),conds,pref_cons,contrasts,dims,paths,coreg_params,1,'_vs-silence'); % clust threshold
compute_clustered_betas(p_val_peak(1),conds,pref_cons,contrasts,dims,paths,coreg_params,0,'_vs-silence'); % peak threshold

compute_clustered_betas(p_val_peak(6),conds,pref_cons,contrasts,dims,paths,coreg_params,0,'_vs-silence'); % peak threshold
compute_clustered_betas(p_val_peak(7),conds,pref_cons,contrasts,dims,paths,coreg_params,0,'_vs-silence'); % peak threshold
compute_clustered_betas(p_val_peak(8),conds,pref_cons,contrasts,dims,paths,coreg_params,0,'_vs-silence'); % peak threshold
compute_clustered_betas(p_val_peak(end),conds,pref_cons,contrasts,dims,paths,coreg_params,0,'_vs-silence'); % peak threshold






% %% Compute prefs from betas
% conds = nconds+1:nconds+nconds-1;
% compute_clustered_betas(conds,pref_cons,contrasts,dims,paths,coreg_params,1); % clust threshold
% compute_clustered_betas(conds,pref_cons,contrasts,dims,paths,coreg_params,0); % peak threshold
