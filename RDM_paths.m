%% Init
cd /hpc/banco/Primavoice_Scripts/Pipeline3
parameters_sub_07 % run a parameter file to add everything needed in the path




% To use with plot_mean_brain_RDMs:
%% Macaque A1
opt.clust_names = {'M-A1-L';'M-A1-R'};
opt.clust_nums = [1 2;1 2]; % LR 19
opt.stim_names = {'HspF';'HspM';'HnsF';'HnsM';'Mcoo';'Mgru';'Magg';'Mscr';'MTtr';'MTph';'MTtw';'MTts';'NVnL';'NVnN';'NVaL';'NVaN'};
% L2 19
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/FSL_mcflirt_dof6_flirt_dof6_R2s_67/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L2/sub-Maga_res-7WM_4CSF_0mvt_euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L2.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Apache/FSL_mcflirt_dof6_flirt_dof6_R2s_64/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L2/sub-Apache_res-7WM_4CSF_0mvt_euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L2.mat';

% L5 19
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/FSL_mcflirt_dof6_flirt_dof6_R2s_67/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L5/sub-Maga_res-7WM_4CSF_0mvt_euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L5.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Apache/FSL_mcflirt_dof6_flirt_dof6_R2s_64/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L5/sub-Apache_res-7WM_4CSF_0mvt_euclidean_brain_RDM_petkov_a1_thres-0p5_clust-size-30_n-vox-19_n-peaks-1_L5.mat';



%% Macaque aTVA
opt.clust_names = {'M-aTVA-L';'M-aTVA-R'};
opt.clust_nums = [2 1;4 7]; % LR 19
% opt.clust_nums = [2 1;15 21]; % LR 19, same L different R, to use with second Apache RDM file
% opt.clust_nums = [2 1;5 21]; % LR 19, different L & R, to use with second Apache RDM file
opt.stim_names = {'HspF';'HspM';'HnsF';'HnsM';'Mcoo';'Mgru';'Magg';'Mscr';'MTtr';'MTph';'MTtw';'MTts';'NVnL';'NVnN';'NVaL';'NVaN'};
% L2 19
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/FSL_mcflirt_dof6_flirt_dof6_R2s_67/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_macaque_vs_all_thres-4p9_clust-size-20_n-vox-19_n-peaks-1_L2/sub-Maga_res-7WM_4CSF_0mvt_euclidean_brain_RDM_macaque_vs_all_thres-4p9_clust-size-20_n-vox-19_n-peaks-1_L2.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Apache/FSL_mcflirt_dof6_flirt_dof6_R2s_64/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_macaque_vs_nonvoice_cut_thres-2_clust-size-30_n-vox-19_n-peaks-0p2_L2/sub-Apache_res-7WM_4CSF_0mvt_euclidean_brain_RDM_macaque_vs_nonvoice_cut_thres-2_clust-size-30_n-vox-19_n-peaks-0p2_L2.mat';
% RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Apache/FSL_mcflirt_dof6_flirt_dof6_R2s_64/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_macaque_vs_all_thres-2_clust-size-30_n-vox-19_n-peaks-0_L2/sub-Apache_res-7WM_4CSF_0mvt_euclidean_brain_RDM_macaque_vs_all_thres-2_clust-size-30_n-vox-19_n-peaks-0_L2.mat';
% L5 19
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/FSL_mcflirt_dof6_flirt_dof6_R2s_67/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_macaque_vs_all_thres-4p9_clust-size-20_n-vox-19_n-peaks-1_L5/sub-Maga_res-7WM_4CSF_0mvt_euclidean_brain_RDM_macaque_vs_all_thres-4p9_clust-size-20_n-vox-19_n-peaks-1_L5.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Apache/FSL_mcflirt_dof6_flirt_dof6_R2s_64/results_7WM_4CSF_0mvt/brain_RDMs/euclidean_brain_RDM_macaque_vs_nonvoice_cut_thres-2_clust-size-30_n-vox-19_n-peaks-0p2_L5/sub-Apache_res-7WM_4CSF_0mvt_euclidean_brain_RDM_macaque_vs_nonvoice_cut_thres-2_clust-size-30_n-vox-19_n-peaks-0p2_L5.mat';

opt.colormap = 'hot';
% opt.stim_names = {'H';'M';'MT';'NV'};


coordsL = [-25 2 -3;-23 1 -10.5];
mean_coordsL = [-24 1.5 -6.75];
coordsR = [24.5 -1.5 -3;20 -3 -8];
mean_coordsR = [22.5 -2.25 -5.5];






%% Human A1 19
opt.clust_names = {'H-A1-L';'H-A1-R'};
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-02/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2/sub-02_res-8WM_9CSF_0mvt_euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-04/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2/sub-04_res-8WM_9CSF_0mvt_euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2.mat';
RDM_files{3} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-05/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2/sub-05_res-8WM_9CSF_0mvt_euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2.mat';
RDM_files{4} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-06/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2/sub-06_res-8WM_9CSF_0mvt_euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2.mat';
RDM_files{5} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-07/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2/sub-07_res-8WM_9CSF_0mvt_euclidean_brain_RDM_heschl_prob_thres-30_clust-size-30_n-vox-19_n-peaks-1_L2.mat';



%% Human TVAs 19
opt.clust_names = {'H-aTVA-L';'H-aTVA-R'};
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-02/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2/sub-02_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-04/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2/sub-04_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2.mat';
RDM_files{3} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-05/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2/sub-05_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2.mat';
RDM_files{4} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-06/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2/sub-06_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2.mat';
RDM_files{5} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-07/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2/sub-07_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-19_n-peaks-0_L2.mat';

opt.clust_nums = [2 7;3 5;7 3;9 5;7 8]; % LR aTVA 19
opt.distances = [4 4;8 15;2 6;8 4;7 9]; % distance from Aglieri (rounded mm)
coordsL = [-59 -3 -3;-56 -2 -2;-64 -4 -3;-54 -6 2;-61 1 -6];
mean_coordsL = [-58.8 -2.8 -2.4];
coordsR = [59 6 -8;55 14 -16;61 9 -5;60 0 -4;60 7 -1];
mean_coordsR = [59 7.2 -6.8];


opt.clust_names = {'H-mTVA-L';'H-mTVA-R'};
opt.clust_nums = [1 6;2 1;9 8;1 11;1 2]; % LR mTVA
opt.distances = [11 9;10 17;10 5;12 8;17 19]; % distance from Aglieri (rounded mm)

opt.clust_names = {'H-pTVA-L';'H-pTVA-R'};
opt.clust_nums = [5 8;6 4;11 1;4 2;3 6]; % LR pTVA
opt.distances = [12 5;5 19;8 6;9 20;18 13]; % distance from Aglieri (rounded mm)

% TVAs combined
opt.clust_names = {'H-aTVA-L';'H-aTVA-R';'H-mTVA-L';'H-mTVA-R';'H-pTVA-L';'H-pTVA-R'};
opt.clust_nums = [2 7 1 6 5 8;3 5 2 1 6 4;7 3 9 8 11 1;9 5 1 11 4 2;7 8 1 2 3 6];

%% FVAs
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-02/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2/sub-02_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2.mat';
RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-04/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2/sub-04_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2.mat';
RDM_files{3} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-05/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-2_clust-size-20_n-vox-19_n-peaks-0_L2/sub-05_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-2_clust-size-20_n-vox-19_n-peaks-0_L2.mat';
RDM_files{4} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-06/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2/sub-06_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-4_clust-size-20_n-vox-19_n-peaks-0_L2.mat';
RDM_files{5} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-07/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-2_clust-size-20_n-vox-19_n-peaks-0_L2/sub-07_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-2_clust-size-20_n-vox-19_n-peaks-0_L2.mat';

opt.clust_names = {'H-aFVA-L';'H-aFVA-R'};
opt.clust_nums = [18 29;11 10;11 13;13 15;17 32]; % LR aFVA

opt.clust_names = {'H-mFVA-L';'H-mFVA-R'};
opt.clust_nums = [13 23;12 15;16 25;21 11;18 25]; % LR mFVA

opt.clust_names = {'H-pFVA-L';'H-pFVA-R'};
opt.clust_nums = [16 32;21 23;29 24;14 16;22 15]; % LR pFVA

% FVAs combined
opt.clust_names = {'H-aFVA-L';'H-aFVA-R';'H-mFVA-L';'H-mFVA-R';'H-pFVA-L';'H-pFVA-R'};
opt.clust_nums = [18 29 13 23 16 32;11 10 12 15 21 23;11 13 16 25 29 24;13 15 21 11 14 16;17 32 18 25 22 15];


%% Human TVAs 57
% RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-02/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2/sub-02_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2.mat';
% RDM_files{2} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-04/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2/sub-04_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2.mat';
% RDM_files{3} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-05/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2/sub-05_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2.mat';
% RDM_files{4} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-06/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2/sub-06_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2.mat';
% RDM_files{5} = '/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-07/spm_realign/results_8WM_9CSF_0mvt/brain_RDMs/euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2/sub-07_res-8WM_9CSF_0mvt_euclidean_brain_RDM_human_vs_all_thres-7_clust-size-300_n-vox-57_n-peaks-0_L2.mat';
% opt.clust_names = {'H-aTVA-L';'H-aTVA-R'};
% opt.clust_nums = [2 3;2 3;2 3;2 3;2 3]; % LR aTVA 57

% % combined
% opt.clust_nums = [2 7 1 6 5 8;3 5 2 1 6 4;7 3 9 8 11 1;9 5 1 11 4 2;7 8 1 2 3 6]; % LR aTVA mTVA pTVA
% opt.distances = [4 4 11 9 12 5;8 15 10 17 5 19;2 6 10 17 8 6;8 4 12 8 9 20;7 9 17 19 18 13]; % distance from Aglieri (rounded mm)
% opt.clust_names = {'Left aTVA';'Left mTVA';'Left pTVA';'Right aTVA';'Right mTVA';'Right pTVA'};


%% Fixed Effect Human
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/second_level_analysis_human/brain_RDMs/fixed_euclidean_brain_RDM_human_vs_all_thres-3p09_clust-size-300_n-vox-297_n-peaks-0_L2/second_level_human_fixed_euclidean_brain_RDM_human_vs_all_thres-3p09_clust-size-300_n-vox-297_n-peaks-0_L2.mat';
opt.clust_names = {'Left aTVA';'Left mTVA';'Left pTVA';'Right aTVA';'Right mTVA';'Right pTVA';'Left aFVA';'Left mFVA';'Left pFVA';'Right aFVA';'Right mFVA';'Right pFVA'};
opt.clust_nums = [7 5 11 4 3 6 30 17 20 21 23 22];
opt.clust_nums = [7 5 8  4 3 6 13 17 20 29 23 22];
opt.coords =     [-64 3 -6;-67 -28 6;-54 -45 7 ;60 6 -6;65 -16 -3;56 -30 4;-34 26 5;-51 11 15;-54 -2 54;57 39 2;49 17 23;53 -1 51];
opt.coords =     [-64 3 -6;-67 -28 6;-66 -43 14;60 6 -6;65 -16 -3;56 -30 4;-54 30 4;-51 11 15;-54 -2 54;50 24 -4;49 17 23;53 -1 51];
aglieri_coords = [-62 -4 0;-66 -28 4;-58 -38 6 ;58 2 -8;58 -20 -2;50 -32 4;-39 27 -3;-48 16 21;-52 -8 48;54 32 0;48 18 24;51 -2 48];

opt.distances = zeros(1,length(opt.clust_names));
for i = 1:length(opt.clust_names)
	opt.distances(i) = norm(aglieri_coords(i,:) - opt.coords(i,:));
end

% opt.stim_names = {'H';'M';'MT';'NV'};
opt.stim_names = {'HspF';'HspM';'HnsF';'HnsM';'Mcoo';'Mgru';'Magg';'Mscr';'MTtr';'MTph';'MTtw';'MTts';'NVnL';'NVnN';'NVaL';'NVaN'};


%% Fixed Effect Macaques
RDM_files{1} = '/hpc/banco/Primavoice_Data_and_Analysis/second_level_analysis_macaque/brain_RDMs/fixed_euclidean_brain_RDM_macaque_vs_all_thres-3p09_clust-size-513_n-vox-513_n-peaks-0_CV-macaques_L2/second_level_macaque_fixed_euclidean_brain_RDM_macaque_vs_all_thres-3p09_clust-size-513_n-vox-513_n-peaks-0_CV-macaques_L2.mat';



plot_mean_brain_RDMs(RDM_files,opt)







%% Save RDMs
opt.dispFigs = 0; % to output RDMs without display figures
[MA1,kMA1,mag_MA1] = plot_mean_brain_RDMs(RDM_files,opt);
[HA1,kHA1,mag_HA1] = plot_mean_brain_RDMs(RDM_files,opt);
[MaTVA,kMaTVA,mag_MaTVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HaTVA,kHaTVA,mag_HaTVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HmTVA,kHmTVA,mag_HmTVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HpTVA,kHpTVA,mag_HpTVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HaFVA,kHaFVA,mag_HaFVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HmFVA,kHmFVA,mag_HmFVA] = plot_mean_brain_RDMs(RDM_files,opt);
[HpFVA,kHpFVA,mag_HpFVA] = plot_mean_brain_RDMs(RDM_files,opt);
H_FFX1 = plot_mean_brain_RDMs(RDM_files,opt);
H_FFX2 = plot_mean_brain_RDMs(RDM_files,opt);

modelRDMs

paths.Mean_RDMs = fullfile(paths.dataset,'Mean_RDMs');
filename = fullfile(paths.Mean_RDMs,'L2_MA1_MaTVA_HA1_HaTVA_Models_19-vox.mat');
% filename = fullfile(paths.Mean_RDMs,'L2_MA1_MaTVA_HA1_HaTVA_Models_57-vox.mat');
% save(filename,'MA1','MaTVA','HA1','HaTVA','kMA1','kMaTVA','kHA1','kHaTVA','mag_MA1','mag_MaTVA','mag_HA1','mag_HaTVA','Models');
filename = fullfile(paths.Mean_RDMs,'L2_All_Voice_Areas_Models_19-vox.mat');
save(filename,'MA1','kMA1','mag_MA1','MaTVA','kMaTVA','mag_MaTVA',...
	'HA1','kHA1','mag_HA1','HaTVA','kHaTVA','mag_HaTVA',...
	'HmTVA','kHmTVA','mag_HmTVA','HpTVA','kHpTVA','mag_HpTVA',...
	'HaFVA','kHaFVA','mag_HaFVA','HmFVA','kHmFVA','mag_HmFVA',...
	'HpFVA','kHpFVA','mag_HpFVA','Models');


load(filename) % load RDMs


%% MDS, Relatedness...
userOptions.projectName = 'project';
userOptions.analysisName = 'analysis';
userOptions.rootPath = fullfile(paths.dataset,'Mean_RDMs');
userOptions.distanceMeasure = 'Spearman';
userOptions.RDMcorrelationType = 'Kendall_taua'; % (which is appropriate whenever categorical models are tested)
userOptions.saveFiguresPS = 0;
userOptions.saveFiguresFig = 0;
userOptions.saveFiguresPDF = 1;
% userOptions.plotpValues = '=';
userOptions.nBootstrap = 100;
userOptions.nRandomisations = 10000;

rsa.MDSRDMs({HA1;MA1;HaTVA;MaTVA;Models},userOptions) % RDM correlations + MDS
rsa.MDSRDMs({MA1;MaTVA;HA1;HaTVA;HmTVA;HpTVA;HaFVA;HmFVA;HpFVA;Models},userOptions) % RDM correlations + MDS
rsa.MDSRDMs({HaTVA;HmTVA;HpTVA;HaFVA;HmFVA;HpFVA},userOptions) % RDM correlations + MDS
rsa.MDSRDMs({H_FFX1(1);H_FFX1(2);H_FFX1(3);H_FFX1(4);H_FFX1(5);H_FFX1(6);H_FFX1(7);H_FFX1(8);H_FFX1(9);H_FFX1(10);H_FFX1(11);H_FFX1(12)},userOptions) % RDM correlations + MDS
rsa.MDSRDMs({H_FFX2},userOptions) % RDM correlations + MDS



userOptions.conditionLabels = {'HspF';'HspM';'HnsF';'HnsM';'Mcoo';'Mgru';'Magg';'Mscr';'MTtr';'MTph';'MTtw';'MTts';'NVnL';'NVnN';'NVaL';'NVaN'};
% userOptions.conditionLabels = {'human';'human';'human';'human';'macaque';'macaque';'macaque';'macaque';'marmoset';'marmoset';'marmoset';'marmoset';'nonvocal';'nonvocal';'nonvocal';'nonvocal'};
ns = 96;
condcols = [50 100 233;233 50 50;233 50 233;233 133 50] ./ 255;

CL = nan(ns,3);
j = 1;
for i = 1:size(condcols,1)
	k = j + (ns / size(condcols,1)) - 1;
	CL(j:k,:) = repmat(condcols(i,:),ns / size(condcols,1),1);
	j = k + 1;
end
userOptions.conditionColours = CL;
userOptions.conditionLabels = repmat({''},ns,1);

rsa.MDSConditions(MaTVA,userOptions)


%% RDMs' relatedness

% Human model vs. brain RDMs
stats_Human=rsa.compareRefRDM2candRDMs(Models(1),{kHA1(1,:);kMA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:)},userOptions);
stats_Human=rsa.compareRefRDM2candRDMs(Models(1),{kHA1(1,:);kMA1(1,:);kMaTVA(1,:);kHaTVA(1,:);kHmTVA(1,:);kHpTVA(1,:);kHaFVA(1,:);kHmFVA(1,:);kHpFVA(1,:);...
	kHA1(2,:);kMA1(2,:);kMaTVA(2,:);kHaTVA(2,:);kHmTVA(2,:);kHpTVA(2,:);kHaFVA(2,:);kHmFVA(2,:);kHpFVA(2,:)},userOptions);

% Macaque model vs. brain RDMs
stats_Macaque=rsa.compareRefRDM2candRDMs(Models(2),{kHA1(1,:);kMA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:)},userOptions);

% Nonvocal model vs. brain RDMs
stats_Nonvocal=rsa.compareRefRDM2candRDMs(Models(3),{kHA1(1,:);kMA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:)},userOptions);


% MA1_L vs. all RDMs
stats_MA1_L=rsa.compareRefRDM2candRDMs(kMA1(1,:),{kHA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);
% MA1_R vs. all RDMs
stats_MA1_R=rsa.compareRefRDM2candRDMs(kMA1(2,:),{kHA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(1,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);


% HA1_L vs. all RDMs
stats_HA1_L=rsa.compareRefRDM2candRDMs(kHA1(1,:),{kMA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);
% HA1_R vs. all RDMs
stats_HA1_R=rsa.compareRefRDM2candRDMs(kHA1(2,:),{kMA1(1,:);kHaTVA(1,:);kMaTVA(1,:);kHA1(1,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);


% HaTVA_L vs. all RDMs
stats_HaTVA_L=rsa.compareRefRDM2candRDMs(kHaTVA(1,:),{kMA1(1,:);kHA1(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);
% HaTVA_R vs. all RDMs
stats_HaTVA_R=rsa.compareRefRDM2candRDMs(kHaTVA(2,:),{kMA1(1,:);kHA1(1,:);kMaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(1,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);


% MaTVA_L vs. all RDMs
stats_MaTVA_L=rsa.compareRefRDM2candRDMs(kMaTVA(1,:),{kMA1(1,:);kHA1(1,:);kHaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(2,:);Models(1);Models(2);Models(3)},userOptions);
% MaTVA_R vs. all RDMs
stats_MaTVA_R=rsa.compareRefRDM2candRDMs(kMaTVA(2,:),{kMA1(1,:);kHA1(1,:);kHaTVA(1,:);kHA1(2,:);kMA1(2,:);kHaTVA(2,:);kMaTVA(1,:);Models(1);Models(2);Models(3)},userOptions);



filename = fullfile(paths.Mean_RDMs,'RDMs_relatedness_stats.mat');
%save(filename,'stats_Human','stats_Macaque','stats_Nonvocal','stats_MA1_L','stats_MA1_R','stats_HA1_L','stats_HA1_R','stats_HaTVA_L','stats_HaTVA_R','stats_MaTVA_L','stats_MaTVA_R')