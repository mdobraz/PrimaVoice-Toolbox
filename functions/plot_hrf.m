% Plot HRF (launch after parameters file)
time=(0:300)/10;
hrf0=fmridesign(time,0,[1 0],[],GLM_params.hrf);
figure
plot(time,squeeze(hrf0.X(:,1,1,1)),'LineWidth',2)
grid on
xlabel('time (seconds)')
ylabel('hrf')
