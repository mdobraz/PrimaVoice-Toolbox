function save_crop_vol(file_to_crop,Y,crop_box,bt_save)

%% Button
bt_save.Enable = 'Off';
bt_save.String = 'Saving...';

pause(0.001) % button update does not work without pause

%% filenames
[pathstr,name,ext] = fileparts(file_to_crop);

if strcmp(ext,'.gz')
   [~,name,ext] = fileparts(name);
end

cropped_nii_file = fullfile(pathstr,[name 'Cropped' ext]);
crop_box_file = fullfile(pathstr,[name 'Cropped.cropbox']);
save(crop_box_file,'crop_box','-ascii'); % save crop_box

%% fslroi parameters
fslroi_params = nan(3,2); % <xmin> <xsize> <ymin> <ysize> <zmin> <zsize>
fslroi_params(1,:) = sort(size(Y,1) - crop_box(1,:));
fslroi_params(2,:) = sort(size(Y,2) - crop_box(2,:));
fslroi_params(3,:) = sort(crop_box(3,:));
fslroi_params(:,2) = fslroi_params(:,2) - fslroi_params(:,1);

%% Find fslroi command
FSL_prefix = 'fsl5.0-';
if exist('paths','var')
    if isfield(paths,'FSL_prefix')
        FSL_prefix = paths.FSL_prefix;
    end
end
[s,~] = system([FSL_prefix 'fslroi']);
if s == 127
    [s,~] = system('fslroi');
    if s == 127
        error('Cannot find FSL fast command \nIf the command needs a prefix (eg. fsl6.0-bet), put this prefix in a ''paths.FSL_prefix'' variable%s','')
    else
        FSL_prefix = '';
    end
end

%% Delete output file if existant (to avoid conflicts between the a compressed and an uncompressed version of the files)
if exist(cropped_nii_file,'file')
    delete(cropped_nii_file);
end

%% Run fslroi
system(sprintf('%sfslroi %s %s %s',FSL_prefix,file_to_crop,cropped_nii_file,sprintf('%i ',fslroi_params')));
system(sprintf('gunzip %s',cropped_nii_file));

%% Display
fprintf('\n\nCropped file saved as:\n')
fprintf('%s\n',cropped_nii_file)
fprintf('Crop box saved as:\n')
fprintf('%s\n\n',crop_box_file)

bt_save.String = 'Saved';
% bt_save.ForegroundColor = 'r';
% bt_save.Enable = 'On';





%% TRASH

% function save_crop_vol(file_to_crop,Ycropped,crop_box,P,bt_save)
% 
% [pathstr,name,ext] = fileparts(file_to_crop);
% 
% if strcmp(ext,'.gz')
%    [~,name,ext] = fileparts(name);
% end
% 
% cropped_file = fullfile(pathstr,[name 'Cropped' ext]);
% P.dim = size(Ycropped);
% P.fname = cropped_file;
% spm_write_vol(P,Ycropped);
% 
% filename = fullfile(pathstr,[name 'Cropped.cropbox']);
% save(filename,'crop_box','-ascii'); % save crop_box
% 
% fprintf('\n\nCropped file saved as:\n')
% fprintf('%s\n',cropped_file)
% fprintf('Crop box saved as:\n')
% fprintf('%s\n\n',filename)
% 
% bt_save.String = 'Saved';
% bt_save.ForegroundColor = 'r';
