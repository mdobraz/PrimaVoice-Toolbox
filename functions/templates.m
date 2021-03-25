% All templates should be in the paths.templates folder
% To add a new template, create a new case in the switch and copy the lines of an existent template and edit the paths


switch paths.template_name
    case 'inia19'
        paths.template = fullfile(paths.templates,'inia19'); % change the last argument to match with the folder containing the template
        paths.template_brain = fullfile(paths.template,'inia19-t1-brain.nii.gz');
        paths.template_brain_mask = fullfile(paths.template,'inia19-brainmask.nii.gz');
        paths.template_head = fullfile(paths.template,'inia19-t1.nii.gz');
        paths.template_tissue.grey = fullfile(paths.template,'inia19-prob_1.nii.gz');
        paths.template_tissue.white = fullfile(paths.template,'inia19-prob_2.nii.gz');
        paths.template_tissue.csf = fullfile(paths.template,'inia19-prob_0.nii.gz');
        
    case 'NMT'
        paths.template = fullfile(paths.templates,'NMT_v1.2');
        paths.template_brain = fullfile(paths.template,'NMT_SS.nii.gz');
        paths.template_brain_mask = fullfile(paths.template,'masks','anatomical_masks','NMT_brainmask.nii.gz');
        paths.template_head = fullfile(paths.template,'NMT.nii.gz');
        paths.template_tissue.grey = fullfile(paths.template,'masks','probabilisitic_segmentation_masks','NMT_segmentation_GM.nii.gz');
        paths.template_tissue.white = fullfile(paths.template,'masks','probabilisitic_segmentation_masks','NMT_segmentation_WM.nii.gz');
        paths.template_tissue.csf = fullfile(paths.template,'masks','probabilisitic_segmentation_masks','NMT_segmentation_CSF.nii.gz');
        
    case 'BSI-NI'
        paths.template = fullfile(paths.templates,'The_Marmoset_MRI_Standard_Brain');
        paths.template_brain = fullfile(paths.template,'Template_T1_brain.nii.gz');
        paths.template_head = fullfile(paths.template,'Template_T1_head.nii.gz');
        paths.template_tissue.grey = fullfile(paths.template,'grey.nii.gz');
        paths.template_tissue.white = fullfile(paths.template,'white.nii.gz');
        paths.template_tissue.csf = fullfile(paths.template,'csf.nii.gz');

    case 'MNI152'
        paths.template = fullfile(paths.templates,'mni_icbm152_nlin_asym_09c');
        paths.template_brain = fullfile(paths.template,'mni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz');
        paths.template_head = fullfile(paths.template,'mni_icbm152_t1_tal_nlin_asym_09c.nii.gz');
        paths.template_tissue.grey = fullfile(paths.template,'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz');
        paths.template_tissue.white = fullfile(paths.template,'mni_icbm152_wm_tal_nlin_asym_09c.nii.gz');
        paths.template_tissue.csf = fullfile(paths.template,'mni_icbm152_csf_tal_nlin_asym_09c.nii.gz');

    otherwise
        error('Template ''%s'' is not available.\nModify your parameters file.',paths.template_name)
end

% paths.STS_STG_mask = fullfile(paths.template,'inia19-STS_STG_dil.nii.gz');
