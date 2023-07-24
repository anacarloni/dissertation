% plot norm L2 of the difference, WF and RFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [WF; RFA1st; RFA2nd];
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
b = bar(1:4,data);%,'bo-','MarkerFaceColor','b','MarkerSize',8)
b(1).FaceColor = [0.4940 0.1840 0.5560];
b(2).FaceColor = [0.4660 0.6740 0.1880];
b(3).FaceColor = [0 0.4470 0.7410];
% hold on
% bar(RFA1st,'r')%,'rs-','MarkerFaceColor','r','MarkerSize',8)
% bar(RFA2nd,'m')%,'m^-','MarkerFaceColor','m','MarkerSize',8)
grid minor
ylim padded
xlim padded
xticks([1 2 3 4])
xticklabels({'$C_{\ell_h}$','$C_{\ell_\alpha}$','$C_{m_h}$','$C_{m_\alpha}$'})
xlabel('Aerodynamic coefficients')
ylabel('L_2 norm')
legend('Walsh function','RFA First form','RFA Second form','location', ...
    'northeast')
set(gca,'FontSize',21)

figname = sprintf('Case%s_L2norm_WFandRFA',caso);
figures_save