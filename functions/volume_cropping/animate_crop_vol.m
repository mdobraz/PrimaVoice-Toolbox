function animate_crop_vol(Ycropped)

get_groot = get(groot);
L = min(get_groot.MonitorPositions(1,3:4)) * 0.5;
b = min(get_groot.MonitorPositions(1,3:4)) * 0.2;

f = figure('Position',[b b L L],...
                  'Toolbar','none',...
                  'MenuBar','none',...
                  'Name','',...
                  'NumberTitle','off',...
                  'resize','off');
              
H_data = axes('Parent',f);
colormap gray

sY = size(Ycropped);

[m,Im] = max(sY(2:3));
r = min(sY(2:3)) / m;
if Im == 1
    H_data.Position = [(1-r)/2 0 r 1];
else
    H_data.Position = [0 (1-r)/2 1 r];
end

for a = 1:size(Ycropped,1)
    imagesc(squeeze(Ycropped(a,:,:)),'Parent',H_data)
    H_data.XTick = [];
    H_data.YTick = [];
    drawnow
end


[m,Im] = max(sY([1 3]));
r = min(sY([1 3])) / m;
if Im == 1
    H_data.Position = [(1-r)/2 0 r 1];
else
    H_data.Position = [0 (1-r)/2 1 r];
end

for b = 1:size(Ycropped,2)
    imagesc(squeeze(Ycropped(:,b,:)),'Parent',H_data)
    H_data.XTick = [];
    H_data.YTick = [];
    drawnow
end



[m,Im] = max(sY([1 2]));
r = min(sY([1 2])) / m;
if Im == 1
    H_data.Position = [(1-r)/2 0 r 1];
else
    H_data.Position = [0 (1-r)/2 1 r];
end

for c = 1:size(Ycropped,3)
    imagesc(squeeze(Ycropped(:,:,c)),'Parent',H_data)
    H_data.XTick = [];
    H_data.YTick = [];
    drawnow
end

close(f)


% for a = 1:size(Ycropped,1)
%     imagesc(squeeze(Ycropped(a,:,:)),'Parent',H_data)
%     H_data.XTick = [];
%     H_data.YTick = [];
%     drawnow
% end
% 
% for b = 1:size(Ycropped,2)
%     imagesc(squeeze(Ycropped(:,b,:)),'Parent',H_data)
%     H_data.XTick = [];
%     H_data.YTick = [];
%     drawnow
% end
% 
% for c = 1:size(Ycropped,3)
%     imagesc(squeeze(Ycropped(:,:,c)),'Parent',H_data)
%     H_data.XTick = [];
%     H_data.YTick = [];
%     drawnow
% end



