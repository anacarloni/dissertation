%% Plot Impulse Responses for All Methods

% plot impulse responses in stairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
stairs(t,y1(:,1,1),'-','LineWidth',3);
hold on
stairs(t2n,y2(:,1,1),'--','LineWidth',3);
stairs(t3n,y3(:,1,1),'-.','LineWidth',3);
stairs(t4n,y4(:,1,1),'-','LineWidth',3);
stairs(tempo,y5(t<tmax,1,1),':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 0 4e-4 8e-4 12e-4])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
set(gcf,'PaperPositionMode','auto')

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
stairs(t,y1(:,1,2),'-','LineWidth',3);
hold on
stairs(t2n,y2(:,1,2),'--','LineWidth',3);
stairs(t3n,y3(:,1,2),'-.','LineWidth',3);
stairs(t4n,y4(:,1,2),'-','LineWidth',3);
stairs(tempo,y5(t<tmax,1,2),':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-2e-4 0 2e-4 4e-4 6e-4])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
set(gcf,'PaperPositionMode','auto')
legend('CFD (DS)',['ERA, r=',num2str(r)],['BPOD, r=',num2str(r)], ...
    ['ERA/OKID, r=',num2str(r)],'Volterra (1st kernel)','location',...
    'northeast')

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
stairs(t,y1(:,2,1),'-','LineWidth',3);
hold on
stairs(t2n,y2(:,2,1),'--','LineWidth',3);
stairs(t3n,y3(:,2,1),'-.','LineWidth',3);
stairs(t4n,y4(:,2,1),'-','LineWidth',3);
stairs(tempo,y5(t<tmax,2,1),':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
set(gcf,'PaperPositionMode','auto')

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
stairs(t,y1(:,2,2),'-','LineWidth',3);
hold on
stairs(t2n,y2(:,2,2),'--','LineWidth',3);
stairs(t3n,y3(:,2,2),'-.','LineWidth',3);
stairs(t4n,y4(:,2,2),'-','LineWidth',3);
stairs(tempo,y5(t<tmax,2,2),':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 -3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
set(gcf,'PaperPositionMode','auto')

% sgtitle('Impulse Responses')
% linkprop(ax, {'XLim','XGrid','YGrid'});
% ax(1).XLim = [0 .1];
% set(gcf,'PaperPositionMode','auto')
% f2.WindowState = 'maximized';
