function realign_plot_Q(Q,titre)

subplot(2,1,1)
plot(Q(:,1:3))
ylabel('mm')
legend('x translation','y translation','z translation')
grid on
title(titre)
subplot(2,1,2)
plot(Q(:,4:6))
xlabel('Volume')
ylabel('degree')
legend('pitch','roll','yaw')
grid on