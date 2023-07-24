%% Plot Impulse Responses for All Methods

% plot impulse responses in curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,1,1),'-','LineWidth',3);
hold on
plot(t2n,y2(:,1,1),'-','LineWidth',3);
plot(t3n,y3(:,1,1),'--','LineWidth',3);
plot(t4n,y4(:,1,1),'-.','LineWidth',3);
plot(tempo,y11_low(t<tmax).*amp_pg_lowAmp,':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 0 4e-4 8e-4 12e-4])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.56,0.3,0.3])
box(ax,'on')
plot(t,y1(:,1,1),'-','LineWidth',3);
hold on
plot(t2n,y2(:,1,1),'-','LineWidth',3);
plot(t3n,y3(:,1,1),'--','LineWidth',3);
plot(t4n,y4(:,1,1),'-.','LineWidth',3);
plot(tempo,y11_low(t<tmax).*amp_pg_lowAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-20e-5 5e-5])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difROM_ImpResp11_r%d',r);
dir = [fileparts(fileparts(pwd)),filesep,'3_output\r=' num2str(r) '\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,1,2),'-','LineWidth',3);
hold on
plot(t2n,y2(:,1,2),'-','LineWidth',3);
plot(t3n,y3(:,1,2),'--','LineWidth',3);
plot(t4n,y4(:,1,2),'-.','LineWidth',3);
plot(tempo,y21_low(t<tmax).*amp_pt_lowAmp,':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-2e-4 0 2e-4 4e-4 6e-4])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
legend('CFD (DS)',['ERA, r=',num2str(r)],['BPOD, r=',num2str(r)], ...
    ['OKID/ERA, r=',num2str(r)],'Volterra','location',...
    'northeast')

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.27,0.56,0.3,0.3])
box(ax,'on')
plot(t,y1(:,1,2),'-','LineWidth',3);
hold on
plot(t2n,y2(:,1,2),'-','LineWidth',3);
plot(t3n,y3(:,1,2),'--','LineWidth',3);
plot(t4n,y4(:,1,2),'-.','LineWidth',3);
plot(tempo,y21_low(t<tmax).*amp_pt_lowAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-20e-5 5e-5])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difROM_ImpResp12_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,2,1),'-','LineWidth',3);
hold on
plot(t2n,y2(:,2,1),'-','LineWidth',3);
plot(t3n,y3(:,2,1),'--','LineWidth',3);
plot(t4n,y4(:,2,1),'-.','LineWidth',3);
plot(tempo,y12_low(t<tmax).*amp_pg_lowAmp,':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.26,0.3,0.3])
box(ax,'on')
plot(t,y1(:,2,1),'-','LineWidth',3);
hold on
plot(t2n,y2(:,2,1),'-','LineWidth',3);
plot(t3n,y3(:,2,1),'--','LineWidth',3);
plot(t4n,y4(:,2,1),'-.','LineWidth',3);
plot(tempo,y12_low(t<tmax).*amp_pg_lowAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-5e-5 10e-5])
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difROM_ImpResp21_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,2,2),'-','LineWidth',3);
hold on
plot(t2n,y2(:,2,2),'-','LineWidth',3);
plot(t3n,y3(:,2,2),'--','LineWidth',3);
plot(t4n,y4(:,2,2),'-.','LineWidth',3);
plot(tempo,y22_low(t<tmax).*amp_pt_lowAmp,':','LineWidth',3);
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 -3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.26,0.3,0.3])
box(ax,'on')
plot(t,y1(:,2,2),'-','LineWidth',3);
hold on
plot(t2n,y2(:,2,2),'-','LineWidth',3);
plot(t3n,y3(:,2,2),'--','LineWidth',3);
plot(t4n,y4(:,2,2),'-.','LineWidth',3);
plot(tempo,y22_low(t<tmax).*amp_pt_lowAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-5e-5 10e-5])
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difROM_ImpResp22_r%d',r);
figures_save

% linkprop(ax, {'XLim','XTick','XGrid','YGrid'});
% ax(1).XTick = [0 .02 .04 .06 .08 .1];
% ax(1).XLim = [0 0.1];
% f3.WindowState = 'maximized';
% figname = 'ImpResp_low_curves';
% hgsave([figname,'.fig'])
% print(figname,'-dpng')
% print(figname,'-depsc')
% print(figname,'-dpdf','-bestfit')
% movefile(strcat(figname,'.png'),strcat(dir,figname,'.png'))
% movefile(strcat(figname,'.fig'),strcat(dir,figname,'.fig'))
% movefile(strcat(figname,'.pdf'),strcat(dir,figname,'.pdf'))
% movefile(strcat(figname,'.eps'),strcat(dir,figname,'.eps'))
