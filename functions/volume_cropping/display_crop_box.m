function crop_lines = display_crop_box(Y,H_crop,H_data,crop_box,crop_lines)

crop_box(4,1) = size(Y,2) + 1 - crop_box(2,1);
crop_box(4,2) = size(Y,2) + 1 - crop_box(2,2);

cb_dim1 = [3 2;3 1;4 1];

cb_dim2 = [1 2 1 1;1 2 2 2;1 1 1 2;2 2 1 2];

d = 1:3;
cb_axes = {'X';'Y';'Z'};
cb_lines = {'Xx1' 'Xx2' 'Xy1' 'Xy2';'Yx1' 'Yx2' 'Yy1' 'Yy2';'Zx1' 'Zx2' 'Zy1' 'Zy2'};

for p = 1:3
    dims = d(d~=p);
    for l = 1:4
        if ~isempty(crop_lines.(cb_lines{p,l}));delete(crop_lines.(cb_lines{p,l}));end
        H_crop.(cb_axes{p}).Color = 'none';
        hold(H_crop.(cb_axes{p}),'on')
        crop_lines.(cb_lines{p,l}) = plot(H_crop.(cb_axes{p}),[crop_box(cb_dim1(p,1),cb_dim2(l,1)) crop_box(cb_dim1(p,1),cb_dim2(l,2))],[crop_box(cb_dim1(p,2),cb_dim2(l,3)) crop_box(cb_dim1(p,2),cb_dim2(l,4))],'w');
        hold(H_crop.(cb_axes{p}),'off')
        H_crop.(cb_axes{p}).XLim = [1 size(Y,dims(2))];
        H_crop.(cb_axes{p}).YLim = [1 size(Y,dims(1))];
        H_crop.(cb_axes{p}).XTick = [];
        H_crop.(cb_axes{p}).YTick = [];
    end
end


imagesc(1,'Parent',H_data)
H_data.XTick = [];
H_data.YTick = [];