%% Plot DS Responses for All ROMs and CFD Response

% plot arbitrary response in curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,yWF4(:,1),'-','LineWidth',3)
hold on
plot(t,yBPOD(:,1),'--','LineWidth',3)
plot(t,yERAOKID(:,1),'-.','LineWidth',3)
plot(t,yERAlowAmp(:,1),':','LineWidth',3)
% xlim([0 .1])
xlim([149.93 150.1])
ylim padded
% yticks([-4e-4 0 4e-4 8e-4 12e-4]);
% xticks([0 .02 .04 .06 .08 .1]);
xlabel('Dimensionless Time')
ylabel('$C_{\ell}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
% legend('CFD (WF4)',['BPOD, r=',num2str(r)],['OKID/ERA, r=',num2str(r)], ...
%     ['ERA, r=',num2str(r)],'location','northeast')
figname = sprintf('DSResp11_r%d',r);
dir = [fileparts(fileparts(pwd)),filesep,'3_output\r=' num2str(r) '\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(t,yWF4(:,2),'-','LineWidth',3)
hold on
plot(t,yBPOD(:,2),'--','LineWidth',3)
plot(t,yERAOKID(:,2),'-.','LineWidth',3)
plot(t,yERAlowAmp(:,2),':','LineWidth',3)
% xlim([0 0.1])
xlim([149.93 150.1])
ylim padded
% yticks([-2e-4 0 2e-4 4e-4 6e-4]);
% xticks([0 .02 .04 .06 .08 .1]);
xlabel('Dimensionless Time')
ylabel('$C_{m}$','Interpreter','latex')
set(gca,'FontSize',21)
grid on
legend('CFD (WF4)',['BPOD, r=',num2str(r)],['OKID/ERA, r=',num2str(r)], ...
    ['ERA, r=',num2str(r)],'location','northeast')
figname = sprintf('DSResp12_r%d',r);
figures_save
