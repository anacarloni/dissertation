% Sobreposição das Partes Reais e Imaginárias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sprintf(' (Case %s)',caso);
figure
plot(kk ,Gr01 ,'k',kk ,Gi01 ,'r','LineWidth',1.5)
xlim([0 3])
ylim([-10 30])
xlabel('\kappa')
ylabel('G_C_l_,_h')
grid on
set(gca,'FontSize',17);

figure
plot(kk ,Gr02 ,'k',kk ,Gi02 ,'r','LineWidth',1.5)
xlim([0 3])
ylim([-10 20])
ylabel('G_C_l_,_\alpha')
xlabel('\kappa')
grid on
set(gca,'FontSize',17);
legend('\Ree (Case 00)','\Imm (Case 00)',strcat('\Ree',aux),strcat('\Imm',aux),'location','northeast')

figure
plot(kk ,Gr03 ,'k',kk ,Gi03 ,'r','LineWidth',1.5)
xlim([0 3])
ylim([-10 4])
ylabel('G_C_m_,_h')
xlabel('\kappa')
grid on
set(gca,'FontSize',17);

figure
plot(kk ,Gr04 ,'k',kk ,Gi04 ,'r','LineWidth',1.5)
xlim([0 3])
ylim([-4.5 1])
ylabel('G_C_m_,_\alpha')
xlabel('\kappa')
grid on
set(gca,'FontSize',17);
