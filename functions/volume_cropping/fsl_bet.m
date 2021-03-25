function fsl_bet(H,file_to_bet,cursor_pos,ed_f,ed_g,pop_opt,paths,force_out_path)

%% GUI things
Or_str = H.String;
H.Enable = 'Off';
H.String = 'Extracting...';
ed_f.Enable = 'Off';
ed_g.Enable = 'Off';
pop_opt.Enable = 'Off';
pause(0.001) % otherwise the above does not work

%% Find fsl bet command
FSL_prefix = 'fsl5.0-';
if exist('paths','var')
    if isfield(paths,'FSL_prefix')
        FSL_prefix = paths.FSL_prefix;
    end
end
[s,~] = system([FSL_prefix 'bet']);
if s == 127
    [s,~] = system('bet');
    if s == 127
        error('Cannot find FSL bet command \nIf the command needs a prefix (eg. fsl6.0-bet), put this prefix in a ''paths.FSL_prefix'' variable%s','')
    else
        FSL_prefix = '';
    end
end

%% Get file names
[pathstr,name,ext] = fileparts(file_to_bet);
if strcmp(ext,'.gz')
   [~,name,ext] = fileparts(name);
end
if exist('force_out_path','var')
    pathstr = force_out_path;
elseif exist('paths','var')
    if isfield(paths,'segmentation_root') % if paths.segmentation_root exists place outputs there
        if ~exist(paths.segmentation_root,'dir');mkdir(paths.segmentation_root);end % create segmentation folder if non-existant
        pathstr = paths.segmentation_root;
    end
end

out_file = fullfile(pathstr,[name '_BET' ext]);
mask_file =  fullfile(pathstr,[name '_BET_mask' ext]);
params_file =  fullfile(pathstr,[name '_BETparameters.json']);

%% Get -f ang -g options from the GUI
fract_int = str2double(ed_f.String);
vert_grad = str2double(ed_g.String);
if (fract_int > 1) || (fract_int < 0)
    error('Fractional intensity must be comprised between 0 and 1')
end
if (vert_grad > 1) || (vert_grad < -1)
    error('Vertical gradient in f must be comprised between -1 and 1')
end

%% Get brain extraction method from the GUI
opt_str = {'';...
           '-R';...
           '-S';...
           '-B';...
           '-Z';...
           '-F';...
           '-A';...
           '-A2'};

opt = opt_str{pop_opt.Value};

%% Open file slector if T2 file is needed (option -A2)
if pop_opt.Value == 8
    uiwait(msgbox('Select T2 file in the next window'));
    if ~isfield(paths,'anat_file')
        paths.anat_file = file_to_bet;
    end
    file_to_load = file_selector('T2*.nii',paths);
    opt = [opt ' ' file_to_load];
end

%% Delete output files if existant (to avoid conflicts between the a compressed and an uncompressed version of the files)
if exist(out_file,'file')
    delete(out_file);
    delete(mask_file);
end
      
%% Execute the commands
system(sprintf('%sbet %s %s -f %f -g %f -m -c %i %i %i %s',FSL_prefix,file_to_bet,out_file,fract_int,vert_grad,cursor_pos(1),cursor_pos(2),cursor_pos(3),opt));
system(sprintf('gunzip %s',out_file));
system(sprintf('gunzip %s',mask_file));

%% Change GUI button string
H.String = 'Extracted';
pause(0.001) % otherwise the above does not work

%% Open fslview
fslview_command = 'fslview';
[s,~] = system(sprintf('%s %s -l Greyscale %s -l Copper &',fslview_command,file_to_bet,out_file));
if s == 127
    error('Could not find %s! \nChange the ''fslview_command'' variable (3 lines above this error), \nto match with the way fslview is called on your system.',fslview_command)
end

%% Save the BET parameters as a json file
jstruct.opt = opt;
jstruct.fract_int = fract_int;
jstruct.vert_grad = vert_grad;
jstruct.cursor_pos = cursor_pos;
spm_jsonwrite(params_file,jstruct)


%% Save figure if the function figurewrite exists
if exist('figurewrite','file')
    view_slice_overlay(file_to_bet,out_file,0)
    figurewrite(fullfile(pathstr,[name '_BET']),[],[],[],1) % Using GLMdenoise function
end

%% Update GUI
H.String = Or_str;
H.Enable = 'On';
ed_f.Enable = 'On';
ed_g.Enable = 'On';
pop_opt.Enable = 'On';




