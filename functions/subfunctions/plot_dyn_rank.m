function plot_dyn_rank(ax,x,y,ranks,max_rank,xl,yl,titre)

plot(x(ranks(1:max_rank)),y(ranks(1:max_rank)),'o','Parent',ax,'MarkerFaceColor','g')
hold(ax,'on')
plot(x(ranks(max_rank+1:end)),y(ranks(max_rank+1:end)),'o','Parent',ax,'MarkerFaceColor','r')
hold(ax,'off')

grid(ax,'on')
xlabel(ax,xl)
ylabel(ax,yl)
title(ax,titre)