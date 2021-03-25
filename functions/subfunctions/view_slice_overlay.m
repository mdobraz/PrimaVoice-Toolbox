function view_slice_overlay(f1,f2,visible,figwidth,alphaval,cm1,cm2,rel_cursor_pos)

if ~exist('visible','var') || visible == 1
    visible = 'on';
else
    visible = 'off';
end

if ~exist('figwidth','var') || isempty(figwidth)
    figwidth = 1500;
end

if ~exist('alphaval','var') || isempty(alphaval)
    alphaval = 0.4;
end

if ~exist('cm1','var') || isempty(cm1)
    cm1 = 'gray';
end

if ~exist('cm2','var') || isempty(cm2)
    cm2 = 'hot';
end

if ~exist('rel_cursor_pos','var') || isempty(rel_cursor_pos)
    rel_cursor_pos = [0.5 0.5 0.5];
end

%% load files
P1 = spm_vol(f1);
Y1 = spm_read_vols(P1);

P2 = spm_vol(f2);
Y2 = spm_read_vols(P2);

%% Check sizes
if ~all(size(Y1) == size(Y2))
    error('Cannot overlay volumes of different sizes')
end

%% Create figure
print_fig = figure('Position',[10 10 figwidth figwidth/3],'Visible',visible,'Color','k');
print_fig.InvertHardcopy = 'off';

%% Create slice axes
borders = 0.001;
sizeX = (1 - (borders * 4)) / 3; % full size of an axis
sizeY = 1 - (borders * 2); % full size of an axis
second_pos = sizeX + (borders * 2);
third_pos = (sizeX * 2) + (borders * 3);
ax_data.X = axes('Parent',print_fig,'Position',[borders borders sizeX sizeY]);
ax_data.Y = axes('Parent',print_fig,'Position',[second_pos borders sizeX sizeY]);
ax_data.Z = axes('Parent',print_fig,'Position',[third_pos borders sizeX sizeY]);

%% Resize slice axes depending on the volume size
sl_dims = {'X';'Y';'Z'};
sY = size(Y1);
d = 1:3;
for dim = 1:3
    dims = d(d~=dim);
    [m,Im] = max(sY(dims));
    r = min(sY(dims)) / m;
    if Im == 1
        pos = [((sizeX-(sizeX*r))/2) 0 r 1];
    else
        pos = [0 ((sizeY-(sizeY*r))/2) 1 r];
    end
    ax_data.(sl_dims{dim}).Position(1:2) = ax_data.(sl_dims{dim}).Position(1:2) + pos(1:2);
    ax_data.(sl_dims{dim}).Position(3:4) = ax_data.(sl_dims{dim}).Position(3:4) .* pos(3:4);
end

%% Create overlay axes
for dim = 1:3
    ax_over.(sl_dims{dim}) = axes('Parent',print_fig,'Position',ax_data.(sl_dims{dim}).Position);
end

%% Display
cursor_pos = [round(size(Y1,1)*rel_cursor_pos(1));round(size(Y1,2)*rel_cursor_pos(2));round(size(Y1,3)*rel_cursor_pos(3))];

for dim = 1:3
    slice = cursor_pos(dim);
    if dim == 1
        alpha = squeeze(Y2(slice,:,:)) > 0;
        ho.(sl_dims{dim}) = imagesc(squeeze(Y2(slice,:,:)),'Parent',ax_over.(sl_dims{dim}));
        imagesc(squeeze(Y1(slice,:,:)),'Parent',ax_data.(sl_dims{dim}))
    elseif dim == 2
        alpha = squeeze(Y2(:,slice,:)) > 0;
        ho.(sl_dims{dim}) = imagesc(squeeze(Y2(:,slice,:)),'Parent',ax_over.(sl_dims{dim}));
        imagesc(squeeze(Y1(:,slice,:)),'Parent',ax_data.(sl_dims{dim}))
    elseif dim == 3
        alpha = squeeze(Y2(:,:,slice)) > 0;
        ho.(sl_dims{dim}) = imagesc(squeeze(Y2(:,:,slice)),'Parent',ax_over.(sl_dims{dim}));
        imagesc(squeeze(Y1(:,:,slice)),'Parent',ax_data.(sl_dims{dim}))
    end
    alpha = alpha .* alphaval;
    ho.(sl_dims{dim}).AlphaData = alpha;
    colormap(ax_data.(sl_dims{dim}),cm1)
    colormap(ax_over.(sl_dims{dim}),cm2)
    
    ax_data.(sl_dims{dim}).XTick = [];
    ax_data.(sl_dims{dim}).YTick = [];
    
    ax_over.(sl_dims{dim}).Color = 'none';
    ax_over.(sl_dims{dim}).XTick = [];
    ax_over.(sl_dims{dim}).YTick = [];
end






