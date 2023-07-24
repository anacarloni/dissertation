% resposta temporal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,reshape(yDS(1,1,:),[1,length(yDS)]),'LineWidth',3)
ylabel('$C_{\ell_h}$','Interpreter','latex')
xlabel('Dimensionless Time')
set(gca,'FontSize',21)
% ylim([-7e-4 13e-4])
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
% yticks([-4e-4 0 4e-4 8e-4 12e-4])
grid on
ylim padded
figname = 'DSy11';
dir = [fileparts(fileparts(pwd)),filesep,'3_output\DS\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,reshape(yDS(1,2,:),[1,length(yDS)]),'LineWidth',3)
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
xlabel('Dimensionless Time')
set(gca,'FontSize',21)
% ylim([-3.5e-4 6.5e-4])
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
% yticks([-2e-4 0 2e-4 4e-4 6e-4])
grid on
ylim padded
figname = 'DSy12';
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,reshape(yDS(2,1,:),[1,length(yDS)]),'LineWidth',3)
ylabel('$C_{m_h}$','Interpreter','latex')
xlabel('Dimensionless Time')
set(gca,'FontSize',21)
% ylim([-3.7e-4 2.2e-4])
xticks([0 .02 .04 .06 .08 .1])
xlim([0 .1])
% yticks([-3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
grid on
ylim padded
figname = 'DSy21';
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,reshape(yDS(2,2,:),[1,length(yDS)]),'LineWidth',3)
ylabel('$C_{m_\alpha}$','Interpreter','latex')
xlabel('Dimensionless Time')
set(gca,'FontSize',21)
% ylim([-4e-4 2.5e-4])
xticks([0 .02 .04 .06 .08 .1])
xlim([0 .1])
% yticks([-4e-4 -3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
grid on
ylim padded
figname = 'DSy22';
figures_save

% linkprop(ax, {'XLim','XTick','XLabel','XGrid','YGrid'});
% ax(1).XTick = [0 .02 .04 .06 .08 .1];
% ax(1).XLim = [0 0.1];
% ax(1).XLabel.String = 'dimensionless time';
% set(gcf,'PaperPositionMode','auto')