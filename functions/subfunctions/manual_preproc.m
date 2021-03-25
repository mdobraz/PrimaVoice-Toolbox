

%% Cropped version of the main anat file
[anat_path,anat_name] = fileparts(anat_file{1});
anat_file_nocrop = fullfile(anat_path,[anat_name '.nocrop']);

if ~exist(anat_file_nocrop,'file')
    anat_file_cropped = spm_BIDS(BIDS,'data','sub',paths.subject,'type',[paths.main_anat_suffix 'Cropped'],'ses',sprintf('%02.0f',paths.main_anat_session),'modality','anat');
    if isempty(anat_file_cropped)
        answer = questdlg(sprintf('No cropped version of the anatomic scan was found.\nDo you want to crop the anatomic scan now?'),'','Yes','No','No & don''t ask again','Yes');
        if strcmp(answer,'Yes')
            paths.anat_file = anat_file{1};
            force_type = ['_' paths.main_anat_suffix '.nii'];
            crop_volume
        elseif strcmp(answer,'No & don''t ask again')
            spm_jsonwrite(anat_file_nocrop,struct('nocrop_file',[anat_name ' will not be cropped. Delete this ''.nocrop'' file if you want to crop it anyway']))
        end
        warning('Please re-run your parameters file')
        return
    else
        paths.anat_file = anat_file_cropped{1};
        paths.main_anat_suffix = [paths.main_anat_suffix 'Cropped'];
    end
else
    paths.anat_file = anat_file{1};
end


%% N4 version of the main anat file
if isfield(paths,'ANTS_path')
    [anat_path,anat_name] = fileparts(paths.anat_file);
    anat_file_noN4 = fullfile(anat_path,[anat_name '.noN4']);
    paths.N4iterations_file = fullfile(anat_path,[anat_name '.N4iterations']);
    
    if ~exist(anat_file_noN4,'file')
        anat_file_N4 = spm_BIDS(BIDS,'data','sub',paths.subject,'type',[paths.main_anat_suffix 'N4'],'ses',sprintf('%02.0f',paths.main_anat_session),'modality','anat');
        if isempty(anat_file_N4)
            %             uiwait(warndlg('No bias field corrected (N4) version of the anatomic scan was found.'))
            %             answer = questdlg('Do you want to run N4 on the anatomic scan now?','','Yes','No','No & don''t ask again','Yes');
            answer = questdlg(sprintf('No bias field corrected (N4) version of the anatomic scan was found.\nDo you want to run N4 on the anatomic scan now?'),'','Yes','No','No & don''t ask again','Yes');
            
            if strcmp(answer,'Yes')
                [paths.anat_file,N4ed,N4iterations] = ANTS_N4(paths.anat_file,paths.ANTS_path);
                paths.main_anat_suffix = [paths.main_anat_suffix 'N4'];
                spm_jsonwrite( paths.N4iterations_file,N4iterations)
            elseif strcmp(answer,'No')
                warning('Please re-run your parameters file')
                return
            end
            if strcmp(answer,'No & don''t ask again') || (exist('N4ed','var') && ~N4ed)
                spm_jsonwrite(anat_file_noN4,struct('noN4_file',[anat_name ' will not be bias field corrected. Delete this ''.noN4'' file if you want to unbias it it anyway']))
            end
        else
            paths.anat_file = anat_file_N4{1};
            paths.main_anat_suffix = [paths.main_anat_suffix 'N4'];
        end
    end
end


%% Denoised version of the main anat file
[anat_path,anat_name] = fileparts(paths.anat_file);
anat_file_noDenoise = fullfile(anat_path,[anat_name '.noDenoise']);

if ~exist(anat_file_noDenoise,'file')
    anat_file_Denoise = spm_BIDS(BIDS,'data','sub',paths.subject,'type',[paths.main_anat_suffix 'Denoised'],'ses',sprintf('%02.0f',paths.main_anat_session),'modality','anat');
    if isempty(anat_file_Denoise)
        %         uiwait(warndlg('No denoised version of the anatomic scan was found.'))
        %         answer = questdlg('Do you want to denoise the anatomic scan now?','','Yes','No','No & don''t ask again','Yes');
        answer = questdlg(sprintf('No denoised version of the anatomic scan was found.\nDo you want to denoise the anatomic scan now?'),'','Yes','No','No & don''t ask again','Yes');
        if strcmp(answer,'Yes')
            [paths.anat_file,denoised] = spm_sanlm(paths.anat_file);
            paths.main_anat_suffix = [paths.main_anat_suffix 'Denoised'];
        elseif strcmp(answer,'No')
            warning('Please re-run your parameters file')
            return
        end
        if strcmp(answer,'No & don''t ask again') || (exist('denoised','var') && ~denoised)
            spm_jsonwrite(anat_file_noDenoise,struct('noDenoise_file',[anat_name ' will not be denoised. Delete this ''.noDenoise'' file if you want to denoise it it anyway']))
        end
    else
        paths.anat_file = anat_file_Denoise{1};
        paths.main_anat_suffix = [paths.main_anat_suffix 'Denoised'];
    end
end