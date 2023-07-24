% plot resposta ao impulso por DS e por primeiro kernel de Volterra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tmax = 0.1;
tempo = t(t<tmax);
plot(tempo,s11(t<tmax),'LineWidth',3)
hold on
plot(tempo,y11_low(t<tmax),'--','LineWidth',3)
plot(tempo,y11_med(t<tmax),'-.','LineWidth',3)
plot(tempo,y11_high(t<tmax),':','LineWidth',3)
xlabel('Dimensionless Time')
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylabel('$C_{\ell_h}$','Interpreter','latex')
ylim padded
set(gca,'FontSize',21)
grid on
% figname = 'ImpResp11_difAmp_Volt';
% dir = [fileparts(fileparts(pwd)),filesep,'3_output\amplitude\'];
% figures_save

figure
plot(tempo,s21(t<tmax),'LineWidth',3)
hold on
plot(tempo,y21_low(t<tmax),'--','LineWidth',3)
plot(tempo,y21_med(t<tmax),'-.','LineWidth',3)
plot(tempo,y21_high(t<tmax),':','LineWidth',3)
xlabel('Dimensionless Time')
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
ylim padded
set(gca,'FontSize',21)
grid on
legend('CFD (DS)','Volterra (Low Amplitude)','Volterra (Medium Amplitude)',...
    'Volterra (High Amplitude)','location','northeast');
% figname = 'ImpResp21_difAmp_Volt';
% figures_save

figure
plot(tempo,s12(t<tmax),'LineWidth',3)
hold on
plot(tempo,y12_low(t<tmax),'--','LineWidth',3)
plot(tempo,y12_med(t<tmax),'-.','LineWidth',3)
plot(tempo,y12_high(t<tmax),':','LineWidth',3)
xlabel('Dimensionless Time')
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylabel('$C_{m_h}$','Interpreter','latex')
ylim padded
set(gca,'FontSize',21)
grid on
% figname = 'ImpResp12_difAmp_Volt';
% figures_save

figure
plot(tempo,s22(t<tmax),'LineWidth',3)
hold on
plot(tempo,y22_low(t<tmax),'--','LineWidth',3)
plot(tempo,y22_med(t<tmax),'-.','LineWidth',3)
plot(tempo,y22_high(t<tmax),':','LineWidth',3)
xlabel('Dimensionless Time')
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylabel('$C_{m_\alpha}$','Interpreter','latex')
ylim padded
set(gca,'FontSize',21)
grid on
% figname = 'ImpResp22_difAmp_Volt';
% figures_save

% linkprop(ax, {'XLim','XTick','XGrid','YGrid'});
% ax(1).XLabel.String = 'Dimensionless Time';
% ax(1).XTick = [0 .02 .04 .06 .08 .1];
% ax(1).XLim = [0 0.1];
% set(gcf,'Position',[700 100 600 550])
% set(gcf,'PaperPositionMode','auto')