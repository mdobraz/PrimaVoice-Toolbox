%% GUI dimensions
borders = 0.08;
slider_size = 0.02;
bt_width = 0.08;
bt_height = slider_size * 1.9;

%% Usefull variables
sl_dims = {'X';'Y';'Z'};
% dim_colors = ['r';'g';'b'];
dim_colors = lines(3);
slice_lines = struct('X1',[],'X2',[],'Y1',[],'Y2',[],'Z1',[],'Z2',[]);

%% Create figure
get_groot = get(groot);
L = min(get_groot.MonitorPositions(1,3:4)) * 0.9;
b = min(get_groot.MonitorPositions(1,3:4)) * 0.1;

if exist('slice_fig','var')
    if isgraphics(slice_fig); delete(slice_fig); end
end
slice_fig = figure('Position',[b b L L],...
                  'Color',[0.9 0.91 0.92],...
                  'Toolbar','none',...
                  'MenuBar','none',...
                  'Name',file_to_load,...
                  'NumberTitle','off',...
                  'resize','off',...
                  'CloseRequestFcn','clear Y Ycropped P; delete(slice_fig)');
%                   'CloseRequestFcn','delete(slice_fig)');

%% Create slice axes
sizeI = (1 - (borders * 3)) / 2; % full size of an axis
second_pos = sizeI + (borders * 2);
ax_data.X = axes('Parent',slice_fig,'Position',[borders second_pos sizeI sizeI]);
ax_data.Y = axes('Parent',slice_fig,'Position',[second_pos second_pos sizeI sizeI]);
ax_data.Z = axes('Parent',slice_fig,'Position',[borders borders sizeI sizeI]);

%% Resize slice axes depending on the volume size
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
end

%% Create slice lines axes
for dim = 1:3
    ax_lines.(sl_dims{dim}) = axes('Parent',slice_fig,'Position',ax_data.(sl_dims{dim}).Position);
end

%% Initial slice positioning
cursor_pos = [round(size(Y,1)/2);round(size(Y,2)/2);round(size(Y,3)/2)];

%% Display initial slices
for dim = 1:3
    slice_lines = display_slice(Y,ax_data,ax_lines,cursor_pos(dim),dim,slice_lines,dim_colors);
end

%% Create navigation sliders
for dim = 1:3
    sl_slice.(sl_dims{dim}) = uicontrol('style','slider',...
              'BackgroundColor',dim_colors(dim,:),...
              'unit','normalized',...
              'Value',cursor_pos(dim),...
              'Min',1,...
              'Position',[ax_data.(sl_dims{dim}).Position(1) (ax_data.(sl_dims{dim}).Position(2) + ax_data.(sl_dims{dim}).Position(4)) ax_data.(sl_dims{dim}).Position(3) slider_size],...
              'Max',size(Y,dim),...
              'SliderStep',[1/(size(Y,dim)-1) 10/(size(Y,dim)-1)],...
              'Callback',sprintf('slice_lines = display_slice(Y,ax_data,ax_lines,sl_slice.%s.Value,%i,slice_lines,dim_colors);cursor_pos = get_cursor_pos(sl_slice);',sl_dims{dim},dim));
end


