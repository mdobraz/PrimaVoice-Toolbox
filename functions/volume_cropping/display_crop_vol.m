function Ycropped = display_crop_vol(Y,crop_box,H_data,posA)

crop_box(1,:) = size(Y,1) + 1 - crop_box(1,:);
crop_box(2,:) = size(Y,2) + 1 - crop_box(2,:);

crop_box = sort(crop_box,2);

Ycropped = Y(crop_box(1,1):crop_box(1,2),crop_box(2,1):crop_box(2,2),crop_box(3,1):crop_box(3,2));


sizeI = posA(4);
sY = size(Ycropped);
d = 1:3;
dim = 1;
dims = d(d~=dim);
[m,Im] = max(sY(dims));
r = min(sY(dims)) / m;
if Im == 1
    pos = [((sizeI-(sizeI*r))/2) 0 r 1];
else
    pos = [0 ((sizeI-(sizeI*r))/2) 1 r];
end
H_data.Position(1:2) = posA(1:2) + pos(1:2);
H_data.Position(3:4) = posA(3:4) .* pos(3:4);



imagesc(squeeze(Ycropped(round(size(Ycropped,1)/2),:,:)),'Parent',H_data)
H_data.XTick = [];
H_data.YTick = [];