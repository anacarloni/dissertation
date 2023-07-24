% Comparação dos Resultados - Função da Frequencia Reduzida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(k,reclh,'bo',k,imclh,'go',kk,real(kCLhfit),'m',kk,imag(kCLhfit),'r', ...
    'linewidth',2,'markersize',13);
xlim([0 3]);
ylabel('$C_{\ell_h}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21)
figname = sprintf('Case%s_IntPolG11',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(k,recla,'bo',k,imcla,'go',kk,real(kCLafit),'m',kk,imag(kCLafit),'r', ...
    'linewidth',2,'markersize',13);
xlim([0 3]);
ylim([-5 15]);
ylabel('$C_{\ell_\alpha}$','interpreter','latex')
xlabel('\kappa')
legend('Real - Data set','Imag - Data set','Real - RFA fit',...
    'Imag - RFA fit','location','northeast')
grid on
set(gca,'FontSize',21)
figname = sprintf('Case%s_IntPolG12',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(k,recmh,'bo',k,imcmh,'go',kk,real(kCMhfit),'m',kk,imag(kCMhfit),'r', ...
    'linewidth',2,'markersize',13);
xlim([0 3]);
ylabel('$C_{m_h}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21)
figname = sprintf('Case%s_IntPolG21',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(k,recma,'bo',k,imcma,'go',kk,real(kCMafit),'m',kk,imag(kCMafit),'r', ...
    'linewidth',2,'markersize',13);
xlim([0 3]);
ylabel('$C_{m_\alpha}$','interpreter','latex')
xlabel('\kappa')
grid on
set(gca,'FontSize',21)
figname = sprintf('Case%s_IntPolG22',caso);
figures_save
