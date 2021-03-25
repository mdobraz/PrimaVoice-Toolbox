function slice = addtext2slice(slice,txtstr,txtval)

fig = figure('Visible','off');
text(0,0.015,txtstr,'FontSize',8,'FontWeight','bold','FontSmoothing','off')
F = getframe(gca,[0 0 size(slice,1) size(slice,2)]);
close(fig)

c = F.cdata(:,:,1);
% [i,j] = find(flipud(c)'==0); % normal way (works with FSL view)
[i,j] = find(rot90(rot90(c'))==0); % mirrored (works with Display of minc tools)
ind = sub2ind(size(slice),i,j);
slice(ind) = txtval;
