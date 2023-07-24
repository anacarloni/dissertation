%% Plot DS Responses for All ROMs and CFD Response

% plot DS responses in curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
yDS11 = reshape(yDS(1,1,:),[1,length(yDS)]);
plot(t,yDS11,'-','LineWidth',3)
hold on
plot(t,yBPOD_pg(:,1),'--','LineWidth',3)
plot(t,yERAOKID_pg(:,1),'-.','LineWidth',3)
plot(t,yERAlowAmp_pg(:,1),':','LineWidth',3)
xlim([0 .1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 0 4e-4 8e-4 12e-4]);
xlabel('Dimensionless Time')
ylabel('$C_{\ell_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
figname = sprintf('DSResp11_r%d',r);
dir = [fileparts(fileparts(pwd)),filesep,'3_output\DS\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
yDS12 = reshape(yDS(1,2,:),[1,length(yDS)]);
plot(t,yDS12,'-','LineWidth',3)
hold on
plot(t,yBPOD_pt(:,1),'--','LineWidth',3)
plot(t,yERAOKID_pt(:,1),'-.','LineWidth',3)
plot(t,yERAlowAmp_pt(:,1),':','LineWidth',3)
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-2e-4 0 2e-4 4e-4 6e-4]);
xlabel('Dimensionless Time')
ylabel('$C_{\ell_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
legend('CFD (DS)',['BPOD, r=',num2str(r)],['OKID/ERA, r=',num2str(r)], ...
    ['ERA, r=',num2str(r)],'location','northeast')
figname = sprintf('DSResp12_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
yDS21 = reshape(yDS(2,1,:),[1,length(yDS)]);
plot(t,yDS21,'-','LineWidth',3)
hold on
plot(t,yBPOD_pg(:,2),'--','LineWidth',3)
plot(t,yERAOKID_pg(:,2),'-.','LineWidth',3)
plot(t,yERAlowAmp_pg(:,2),':','LineWidth',3)
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_h}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
figname = sprintf('DSResp21_r%d',r);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
yDS22 = reshape(yDS(2,2,:),[1,length(yDS)]);
plot(t,yDS22,'-','LineWidth',3)
hold on
plot(t,yBPOD_pt(:,2),'--','LineWidth',3)
plot(t,yERAOKID_pt(:,2),'-.','LineWidth',3)
plot(t,yERAlowAmp_pt(:,2),':','LineWidth',3)
xlim([0 0.1])
xticks([0 .02 .04 .06 .08 .1])
ylim padded
yticks([-4e-4 -3e-4 -2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Dimensionless Time')
ylabel('$C_{m_\alpha}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
figname = sprintf('DSResp22_r%d',r);
figures_save