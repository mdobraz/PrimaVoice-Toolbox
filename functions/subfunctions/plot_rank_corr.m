function [selected_sessions,selected_runs,classement,max_rank,selected] = plot_rank_corr(X,Y,xl,yl,selected)

x = X.classement;
y = Y.classement;
mx = X.max_rank;
my = Y.max_rank;

col = lines(6);

%% plot
figure
plot(x((x<=mx) & (y<=my)),y((x<=mx) & (y<=my)),'ko','MarkerFaceColor',col(5,:)) % in both x & y, green
hold on
plot(x((x>mx) & (y<=my)),y((x>mx) & (y<=my)),'ko','MarkerFaceColor',col(6,:)) % in y & not x, light blue
plot(x((x<=mx) & (y>my)),y((x<=mx) & (y>my)),'ko','MarkerFaceColor',col(3,:)) % in x & not y, light orange
plot(x((x>mx) & (y>my)),y((x>mx) & (y>my)),'ko','MarkerFaceColor',col(2,:)) % in not in x & y, red


grid on
xlabel(xl)
ylabel(yl)
[r,pval] = corr(x,y,'type','spearman');
b = robustfit(x,y);
plot(x,b(1)+b(2)*x,'k','LineWidth',2)
if exist('selected','var')
    plot(x(selected),y(selected),'k.')
    [rS,pvalS] = corr(x(selected),y(selected),'type','spearman');
    b = robustfit(x(selected),y(selected));
    plot(x(selected),b(1)+b(2)*x(selected),'r')
    title([sprintf('Spearman rho = %1.2f ',r) '(\color{red}' sprintf('%1.2f',rS) '\color{black})' sprintf(', pval = %1.4f ',pval) '(\color{red}' sprintf('%1.4f',pvalS) '\color{black})'])
%         \color{black}), pval = %1.4f (\color{red}%1.4f\color{black})',r,rS,pval,pvalS))
else
    title(sprintf('Spearman rho = %1.2f, pval = %1.4f',r,pval))
end

plot([mx+0.5 mx+0.5],[0 my+0.5],'--k')
plot([0 mx+0.5],[my+0.5 my+0.5],'--k')


%% selected_sessions
selected_sessions = X.selected_sessions(ismember(X.selected_sessions,Y.selected_sessions));

selected_runs = cell(length(selected_sessions),1);
if ~isempty(selected_sessions)
    for s = 1:length(selected_sessions)
        sess = selected_sessions(s);
        runsX = X.selected_runs{X.selected_sessions == sess};
        runsY = Y.selected_runs{Y.selected_sessions == sess};
        selected_runs{s} = runsX(ismember(runsX,runsY));
    end
else
    error('No session in common!!!')
end


%% Ranking & max_rank
selected = (x<=mx) & (y<=my);
max_rank = sum(selected);
fprintf('Number of selected runs: %i / %i\n',max_rank,length(x))

classement = mean([x y],2);
[~,jnd] = sort(classement,'ascend');
[~,classement] = sort(jnd,'ascend');










