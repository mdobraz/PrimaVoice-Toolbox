SR_ME = struct_sess_run(BIDS,'all','all',paths.subject,[],'epi','anat');

n_runs = 0;
for s = 1:numel(SR_ME)
    session = SR_ME(s).session;
    runs = SR_ME(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        n_runs = n_runs + 1;
        T2s_file = fullfile(paths.ME,sprintf('sub-%s_ses-%02.0f_run-%02.0f_T2s.mat',paths.subject,session,run));
		T2s = load(T2s_file);
		SR_ME(s).mean_R(r) = T2s.mean_R;
		SR_ME(s).time{r} = datevec(T2s.scan_time);
    end
end

opts = optimset('MaxFunEvals',100000,'MaxIter',100000);

col = lines(n_runs);

figure
hold on
xlabel('Time (minutes)')
ylabel('Relaxation rate: R_2* (s-1)')
t_shift = 50;
start_time = 0;
xt = nan(n_runs,1);
xtl = nan(n_runs,1);
i1 = 1;
mid_sess_R2 = nan(numel(SR_ME),1);
B = nan(numel(SR_ME),2);
ref_time = nan(1,6);
n_sess = 0;
for s = 1:numel(SR_ME)
    session = SR_ME(s).session;
    runs = SR_ME(s).runs;
    n_sess = n_sess + 1;
    ref_time = SR_ME(s).time{1};
    times = zeros(length(runs),1);
    for r = 1:length(runs)
    	times(r) = round(etime(SR_ME(s).time{r},ref_time) / 60);
    end
    times_plot = times + start_time;
    i2 = i1 + length(runs) - 1;
    xt(i1:i2) = times_plot;
    xtl(i1:i2) = times;
    plot(times_plot,SR_ME(s).mean_R,'ko','MarkerFaceColor',col(n_sess,:))
    yt = max(SR_ME(s).mean_R) + (std(SR_ME(s).mean_R) / 2);
    text(times_plot(1),yt,sprintf('ses. %i',session))
    % fit monoexponential decay
    if length(runs) > 1
        y = @(b,x) b(1).*exp(-b(2).*x);
        x = times'; % acquisition times
        yx = SR_ME(s).mean_R; % R2*
        OLS = @(b) sum((y(b,x)-yx).^2);
        B(s,:) = fminsearch(OLS,[yx(1) 0.001],opts);
        plot(times_plot,y(B(s,:),x))
    end
%     % find mid-session R2*
%     mid_sess_R2(s) = y(B(s,:),times(end) / 2);
%     plot([times_plot(1) times_plot(end)],[mid_sess_R2(s) mid_sess_R2(s)],'k')
    start_time = start_time + times(end) + t_shift;
    i1 = i2 + 1;
end

set(gca,'XTick',xt)
set(gca,'XTickLabel',num2str(xtl))
grid on





