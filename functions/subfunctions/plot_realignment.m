%%%%%% Run 'set_parameters' beforehand

%% plot
for s = 1:numel(SR)
    session = SR(s).session;
    runs = SR(s).runs;
    for r = 1:length(runs)
        run = runs(r);
        ref_file = sprintf('%s/ref%s_session%i_bold%i.mat',paths.realign,paths.subject,session,run);
        rp_file = sprintf('%s/rp_%s_session%i_bold%i.txt',paths.rp_files,paths.subject,session,run);
        load(ref_file);
        Q = load(rp_file);
        Q(:,4:6) = rad2deg(Q(:,4:6)); % radians to degrees
        titre = sprintf('Session %i, run %i. %i / %i volumes selected',session,run,sum(selected_vols),length(selected_vols));
        
        figure('units','normalized','outerposition',[0 0 1 1]) % fullscreen figure
        realign_plot_Q(Q,titre)
        hold on
        plot([ref_vol ref_vol],[min(min(Q(:,4:6))) max(max(Q(:,4:6)))])
        text(ref_vol+1,(max(max(Q(:,4:6))) / 2),sprintf('representative\nvolume'))
        plot(~selected_vols .* max(max(Q(:,4:6))),'LineWidth',2)
        hold off
        
        pause
        close
    end
end
