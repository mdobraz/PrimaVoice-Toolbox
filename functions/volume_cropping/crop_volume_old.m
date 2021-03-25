%% Choose volume to crop
startdir = '';
if exist('paths','var')
    if isfield(paths,'anat_file')
        startdir = fileparts(paths.anat_file);
    end
end
[file_to_crop,path_ftl] = uigetfile(fullfile(startdir,'*.nii*'),'Select file to crop');

if file_to_crop == 0
    error('Please select one file to crop')
end

file_to_crop = [path_ftl file_to_crop];

%% Load files
P = spm_vol(file_to_crop);
Y = spm_read_vols(P);
Ycropped = Y;

%% GUI dimensions
borders = 0.08;
slider_size = 0.02;
bt_width = 0.08;
bt_height = slider_size * 1.9;

%% Usefull variables
sl_dims = {'X';'Y';'Z'};
dim_colors = ['r';'g';'b'];
slice_lines = struct('X1',[],'X2',[],'Y1',[],'Y2',[],'Z1',[],'Z2',[]);
crop_lines = struct('Xx1',[],'Xx2',[],'Xy1',[],'Xy2',[],'Yx1',[],'Yx2',[],'Yy1',[],'Yy2',[],'Zx1',[],'Zx2',[],'Zy1',[],'Zy2',[]);

%% Create figure
get_groot = get(groot);
L = min(get_groot.MonitorPositions(1,3:4)) * 0.9;
b = min(get_groot.MonitorPositions(1,3:4)) * 0.05;
crop_fig = figure('Position',[b b L L],...
                  'Color',[0.91 0.95 0.97],...
                  'Toolbar','none',...
                  'MenuBar','none',...
                  'Name',file_to_crop,...
                  'NumberTitle','off',...
                  'CloseRequestFcn','clear Y Ycropped P; delete(crop_fig)',...
                  'resize','off');

%% Create axes
sizeI = (1 - (borders * 3)) / 2; % full size of an axis
second_pos = sizeI + (borders * 2);
ax_data.X = axes('Parent',crop_fig,'Position',[borders second_pos sizeI sizeI]);
ax_data.Y = axes('Parent',crop_fig,'Position',[second_pos second_pos sizeI sizeI]);
ax_data.Z = axes('Parent',crop_fig,'Position',[borders borders sizeI sizeI]);
posA = [second_pos borders sizeI sizeI];
ax_data.A = axes('Parent',crop_fig,'Position',posA);

for dim = 1:3
    ax_lines.(sl_dims{dim}) = axes('Parent',crop_fig,'Position',ax_data.(sl_dims{dim}).Position);
    ax_crop.(sl_dims{dim}) = axes('Parent',crop_fig,'Position',ax_data.(sl_dims{dim}).Position);
end

%% Resize axes depending on the volume size
sY = size(Y);
d = 1:3;
for dim = 1:3
    dims = d(d~=dim);
    [m,Im] = max(sY(dims));
    r = min(sY(dims)) / m;
    if Im == 1
        pos = [((sizeI-(sizeI*r))/2) 0 r 1];
    else
        pos = [0 ((sizeI-(sizeI*r))/2) 1 r];
    end
    ax_data.(sl_dims{dim}).Position(1:2) = ax_data.(sl_dims{dim}).Position(1:2) + pos(1:2);
    ax_data.(sl_dims{dim}).Position(3:4) = ax_data.(sl_dims{dim}).Position(3:4) .* pos(3:4);
    ax_lines.(sl_dims{dim}).Position(1:2) = ax_lines.(sl_dims{dim}).Position(1:2) + pos(1:2);
    ax_lines.(sl_dims{dim}).Position(3:4) = ax_lines.(sl_dims{dim}).Position(3:4) .* pos(3:4);
    ax_crop.(sl_dims{dim}).Position(1:2) = ax_crop.(sl_dims{dim}).Position(1:2) + pos(1:2);
    ax_crop.(sl_dims{dim}).Position(3:4) = ax_crop.(sl_dims{dim}).Position(3:4) .* pos(3:4);
    
    if dim == 1
        ax_data.A.Position(1:2) = ax_data.A.Position(1:2) + pos(1:2);
        ax_data.A.Position(3:4) = ax_data.A.Position(3:4) .* pos(3:4);
    end
end

%% Initial crop box & slice positioning
crop_box = [floor((size(Y,1)+1)/4) floor((size(Y,1)+1)/4)*3;...
            floor((size(Y,2)+1)/4) floor((size(Y,2)+1)/4)*3;...
            floor((size(Y,3)+1)/4) floor((size(Y,3)+1)/4)*3];

slice_dim = [round(size(Y,1))/2;round(size(Y,2))/2;round(size(Y,3))/2];

%% Display initial slices
for dim = 1:3
    slice_lines = display_slice(Y,ax_data,ax_lines,slice_dim(dim),dim,slice_lines,dim_colors);
end

%% Create navigation sliders
for dim = 1:3
    sl_slice.(sl_dims{dim}) = uicontrol('style','slider',...
              'BackgroundColor',dim_colors(dim),...
              'unit','normalized',...
              'Value',slice_dim(dim),...
              'Min',1,...
              'Position',[ax_data.(sl_dims{dim}).Position(1) (ax_data.(sl_dims{dim}).Position(2) + ax_data.(sl_dims{dim}).Position(4)) ax_data.(sl_dims{dim}).Position(3) slider_size],...
              'Max',size(Y,dim),...
              'SliderStep',[1/(size(Y,dim)-1) 10/(size(Y,dim)-1)],...
              'Callback',sprintf('slice_lines = display_slice(Y,ax_data,ax_lines,sl_slice.%s.Value,%i,slice_lines,dim_colors);',sl_dims{dim},dim));
end


%% Display initial crop box
crop_lines = display_crop_box(Y,ax_crop,ax_data.A,crop_box,crop_lines);
             
%% Create cropping sliders
sliders = {'Xy1' 'Xy2' 'Xx1' 'Xx2';'Yy1' 'Yy2' 'Yx1' 'Yx2';'Zy1' 'Zy2' 'Zx1' 'Zx2'};
corresp_sliders = {'Zx2' 'Zx1' 'Yx1' 'Yx2';'Zy1' 'Zy2' 'Xx1' 'Xx2';'Yy1' 'Yy2' 'Xy2' 'Xy1'};
d = 1:3;
f = [1 3];
for dim = 1:3
    dims = d(d~=dim);
    l = 0;
    for i = 1:2 % y or x axis
%         nd = dims(i);
        for j = 1:2 % first or second slider
            l = l + 1;
            if (dim == 1 && i == 1) || (dim == 3 && i == 2)
                calbck = sprintf('sl_crop.%s.Value = size(Y,2) + 1 - sl_crop.%s.Value;',corresp_sliders{dim,l},sliders{dim,l});
            else
                calbck = sprintf('sl_crop.%s.Value = sl_crop.%s.Value;',corresp_sliders{dim,l},sliders{dim,l});
            end
            
            if (dim == 3 && i == 2)
                calbck = [calbck sprintf('crop_box(%i,%i) = round(size(Y,2) + 1 - sl_crop.%s.Value);',dims(i),mod(j,2)+1,sliders{dim,l})];
            else
                calbck = [calbck sprintf('crop_box(%i,%i) = round(sl_crop.%s.Value);',dims(i),j,sliders{dim,l})];
            end
                
            calbck = [calbck 'crop_lines = display_crop_box(Y,ax_crop,ax_data.A,crop_box,crop_lines);;bt_anim.Enable = ''Off'';bt_save.Enable = ''Off''; bt_save.String = ''Save''; bt_save.ForegroundColor = [0 0 0];'];
            
            val = floor((size(Y,dims(i))+1)*f(j)/4);
            
            pos = nan(1,4);
            if i == 1 % y axis
                pos(1) = ax_data.(sl_dims{dim}).Position(1) - (slider_size * j);
                pos(2) = ax_data.(sl_dims{dim}).Position(2);
                pos(3) = slider_size;
                pos(4) = ax_data.(sl_dims{dim}).Position(4);
            else % x axis
                pos(1) = ax_data.(sl_dims{dim}).Position(1);
                pos(2) = ax_data.(sl_dims{dim}).Position(2) - (slider_size * j);
                pos(3) = ax_data.(sl_dims{dim}).Position(3);
                pos(4) = slider_size;
            end
            
            sl_crop.(sliders{dim,l}) = uicontrol('style','slider',...
                   'unit','normalized',...
                   'Value',val,...
                   'Min',1,...
                   'Max',size(Y,dims(i)),...
                   'Position',pos,...
                   'SliderStep',[1/(size(Y,dims(i))-1) 10/(size(Y,dims(i))-1)],...
                   'Callback',calbck);
        end
    end
end


%% Create buttons
bt_pos1 = second_pos + (sizeI / 6) - (bt_width / 2);
bt_pos2 = second_pos + (sizeI / 2) - (bt_width / 2);
bt_pos3 = second_pos + (5 * sizeI / 6) - (bt_width / 2);
bt_posy = borders - slider_size*2;

bt_crop = uicontrol('Style','pushbutton',...
                 'units','normalized',...
                 'FontWeight','bold',...
                 'FontSize',(get_groot.MonitorPositions(1,3) / 256) * 1.2,...
                 'String','Crop',...
                 'Enable','on',...
                 'Position',[bt_pos1 bt_posy bt_width bt_height],...
                 'Callback','Ycropped = display_crop_vol(Y,crop_box,ax_data.A,posA);bt_anim.Enable = ''On'';bt_save.Enable = ''On'';');
             
bt_anim = uicontrol('Style','pushbutton',...
                 'units','normalized',...
                 'FontWeight','bold',...
                 'FontSize',(get_groot.MonitorPositions(1,3) / 256) * 1.2,...
                 'String','Animate',...
                 'Enable','off',...
                 'Position',[bt_pos2 bt_posy bt_width bt_height],...
                 'Callback','animate_crop_vol(Ycropped)');
             
bt_save = uicontrol('Style','pushbutton',...
                 'units','normalized',...
                 'FontWeight','bold',...
                 'FontSize',(get_groot.MonitorPositions(1,3) / 256) * 1.2,...
                 'String','Save',...
                 'Enable','off',...
                 'Position',[bt_pos3 bt_posy bt_width bt_height],...
                 'Callback','save_crop_vol(file_to_crop,Ycropped,P,bt_save)'); 








