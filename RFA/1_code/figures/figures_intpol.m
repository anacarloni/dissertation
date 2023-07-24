% Sobreposição das Partes Reais e Imaginárias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sprintf(' (Case %s)',caso);
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(kk_c00,Gr01_c00,'k',kk_c00,Gi01_c00,'r','LineWidth',3)
hold on
plot(k,Gr11,'b',k,Gi11,'g','LineWidth',3)
xlim([0 3])
ylim([-10 30])
xlabel('\kappa')
ylabel('$C_{\ell_h}$','interpreter','latex')
grid on
set(gca,'FontSize',21);
dir = [fileparts(pwd),filesep,'3_output/output_spectral/'];
figname = sprintf('Case%s_TF_G11',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(kk_c00,Gr02_c00,'k',kk_c00,Gi02_c00,'r','LineWidth',3)
hold on
plot(k,Gr12,'b',k,Gi12,'g','LineWidth',3)
xlim([0 3])
ylim([-10 20])
ylabel('$C_{\ell_\alpha}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21);
legend('\Ree (Case 00)','\Imm (Case 00)',strcat('\Ree',aux),strcat('\Imm',aux),'location','northeast')
figname = sprintf('Case%s_TF_G12',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(kk_c00,Gr03_c00,'k',kk_c00,Gi03_c00,'r','LineWidth',3)
hold on
plot(k,Gr21,'b',k,Gi21,'g','LineWidth',3)
xlim([0 3])
ylim([-10 4])
ylabel('$C_{m_h}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21);
% legend('\Ree (Case 00)','\Imm (Case 00)',strcat('\Ree',aux),strcat('\Imm',aux),'location','southwest')
figname = sprintf('Case%s_TF_G21',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(kk_c00,Gr04_c00,'k',kk_c00,Gi04_c00,'r','LineWidth',3)
hold on
plot(k,Gr22,'b',k,Gi22,'g','LineWidth',3)
xlim([0 3])
ylim([-4.5 1])
ylabel('$C_{m_\alpha}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21);
figname = sprintf('Case%s_TF_G22',caso);
figures_save
