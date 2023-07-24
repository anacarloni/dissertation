%% plot ERA Response for Different Amplitudes

% plot impulse responses in plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,1,1)./amp_pg_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,1,1)./amp_pg_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,1,1)./amp_pg_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,1,1)./amp_pg_highAmp,':','LineWidth',3);
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.56,0.3,0.3])
box(ax,'on')
plot(t,y1(:,1,1)./amp_pg_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,1,1)./amp_pg_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,1,1)./amp_pg_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,1,1)./amp_pg_highAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-250 50])
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difAmp_ImpResp11_ERA_r%d',r);
dir = [fileparts(fileparts(pwd)),filesep,'3_output\r=' num2str(r) '\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,1,2)./amp_pt_lowAmp,'-','LineWidth',3);
hold on
plot(t2n,y2(:,1,2)./amp_pt_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,1,2)./amp_pt_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,1,2)./amp_pt_highAmp,':','LineWidth',3);
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
xlabel('Dimensionless Time')
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
legend('CFD (DS)','ERA (Low Amplitude)','ERA (Medium Amplitude)',...
    'ERA (High Amplitude)','location','northeast')
figname = sprintf('difAmp_ImpResp12_ERA_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,2,1)./amp_pg_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,2,1)./amp_pg_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,2,1)./amp_pg_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,2,1)./amp_pg_highAmp,':','LineWidth',3);
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.26,0.3,0.3])
box(ax,'on')
plot(t,y1(:,2,1)./amp_pg_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,2,1)./amp_pg_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,2,1)./amp_pg_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,2,1)./amp_pg_highAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-50 100])
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difAmp_ImpResp21_ERA_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,y1(:,2,2)./amp_pt_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,2,2)./amp_pt_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,2,2)./amp_pt_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,2,2)./amp_pt_highAmp,':','LineWidth',3);
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
set(gcf,'PaperPositionMode','auto')


% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.26,0.3,0.3])
box(ax,'on')
plot(t,y1(:,2,2)./amp_pt_lowAmp,'LineWidth',3);
hold on
plot(t2n,y2(:,2,2)./amp_pt_lowAmp,'--','LineWidth',3);
plot(t3n,y3(:,2,2)./amp_pt_medAmp,'-.','LineWidth',3);
plot(t4n,y4(:,2,2)./amp_pt_highAmp,':','LineWidth',3);
set(ax,'xlim',[.01 .04],'ylim',[-.5 1])
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
xticks([.02 .04])
set(gca,'FontSize',14)
grid minor

figname = sprintf('difAmp_ImpResp22_ERA_r%d',r);
figures_save
