function slice_lines = display_slice(Y,H_data,H_lines,slice,dim,slice_lines,dim_colors)

slice = round(slice);

sl_axes = {'X';'Y';'Z'};
sl_lines = {'X1' 'X2';'Y1' 'Y2';'Z1' 'Z2'};

if dim == 1
    imagesc(squeeze(Y(slice,:,:)),'Parent',H_data.(sl_axes{dim}))
elseif dim == 2
    imagesc(squeeze(Y(:,slice,:)),'Parent',H_data.(sl_axes{dim}))
elseif dim == 3
    imagesc(squeeze(Y(:,:,slice)),'Parent',H_data.(sl_axes{dim}))
end
colormap gray
H_data.(sl_axes{dim}).XTick = [];
H_data.(sl_axes{dim}).YTick = [];
% H_lines.(sl_axes{dim}).Color = 'none';
d = 1:3;
dims = d(d~=dim);

for i = 1:length(dims)
    p = dims(i);
    nodims = d(d~=p);
    l = find(nodims==dim);
    % display slices in other patchlines
    H_lines.(sl_axes{p}).Color = 'none';
    hold(H_lines.(sl_axes{p}),'on')
    if ~isempty(slice_lines.(sl_lines{p,l}));delete(slice_lines.(sl_lines{p,l}));end
    
    disp_slice = slice;
    if l == 1
        if dim == 1; disp_slice = size(Y,dim) - slice + 1;end
        if dim == 2; disp_slice = size(Y,dim) - slice + 1;end
        slice_lines.(sl_lines{p,l}) = plot(H_lines.(sl_axes{p}),[1 size(Y,nodims(2))],[disp_slice disp_slice],'Color',dim_colors(dim,:),'LineWidth',2);
    elseif l == 2
        slice_lines.(sl_lines{p,l}) = plot(H_lines.(sl_axes{p}),[disp_slice disp_slice],[1 size(Y,nodims(1))],'Color',dim_colors(dim,:),'LineWidth',2);
    end
    slice_lines.(sl_lines{p,l}).Color(4) = 0.5; % set alpha value for transparency
    hold(H_lines.(sl_axes{p}),'off')
    H_lines.(sl_axes{p}).XLim = [1 size(Y,nodims(2))];
    H_lines.(sl_axes{p}).YLim = [1 size(Y,nodims(1))];
    H_lines.(sl_axes{p}).XTick = [];
    H_lines.(sl_axes{p}).YTick = [];
end





