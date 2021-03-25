%% Choose volume to crop
if ~exist('force_type','var')
    file_type = '.nii';
else
    file_type = force_type;
    clear force_type
end

if exist('force_in_file','var')
    file_to_load = file_selector(file_type,force_in_file);
    clear force_in_file
elseif exist('paths','var')
    file_to_load = file_selector(file_type,paths);
else
    file_to_load = file_selector(file_type);
end


%% Load volume
P = spm_vol(file_to_load);
Y = spm_read_vols(P);

%% Allocate Ycropped
Ycropped = Y;

%% Display slice viewer
slice_viewer

%% Usefull variables
crop_lines = struct('Xx1',[],'Xx2',[],'Xy1',[],'Xy2',[],'Yx1',[],'Yx2',[],'Yy1',[],'Yy2',[],'Zx1',[],'Zx2',[],'Zy1',[],'Zy2',[]);

%% Create axes
posA = [second_pos borders sizeI sizeI];
ax_data.A = axes('Parent',slice_fig,'Position',posA);

for dim = 1:3
    ax_crop.(sl_dims{dim}) = axes('Parent',slice_fig,'Position',ax_data.(sl_dims{dim}).Position);
end


%% Initial crop box
crop_box = [floor((size(Y,1)+1)/4) floor((size(Y,1)+1)/4)*3;...
            floor((size(Y,2)+1)/4) floor((size(Y,2)+1)/4)*3;...
            floor((size(Y,3)+1)/4) floor((size(Y,3)+1)/4)*3];


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
bt_pos4 = 0.5 - bt_width;
bt_pos5 = 1 - (bt_height / 2);
bt_posy = borders - slider_size*2;

bt_crop = uicontrol('Style','pushbutton',...
                 'units','normalized',...
                 'FontWeight','bold',...
                 'FontSize',(get_groot.MonitorPositions(1,3) / 256) * 1.2,...
                 'String','Crop',...
                 'Enable','on',...
                 'Position',[bt_pos1 bt_posy bt_width bt_height],...
                 'Callback','Ycropped = display_crop_vol(Y,crop_box,ax_data.A,posA);bt_anim.Enable = ''On'';bt_save.Enable = ''On'';;bt_save.String = ''Save'';');
             
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
                 'Callback','save_crop_vol(file_to_load,Y,crop_box,bt_save)'); 
%                  'Callback','save_crop_vol(file_to_load,Ycropped,crop_box,P,bt_save)'); 


bt_loadCB = uicontrol('Style','pushbutton',...
                 'units','normalized',...
                 'FontWeight','bold',...
                 'FontSize',(get_groot.MonitorPositions(1,3) / 256) * 1.2,...
                 'String','Load Crop Box',...
                 'Position',[bt_pos4 bt_pos5 bt_width*2 bt_height/2],...
                 'Callback','cropbox_selector;crop_box = load(CB_to_load);crop_lines = display_crop_box(Y,ax_crop,ax_data.A,crop_box,crop_lines);'); 





