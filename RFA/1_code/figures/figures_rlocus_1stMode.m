% Plotanto os Local das Raízes do Primeiro Modo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
Ax(1) = axes(f2);
L1 = plot(sigma_ref,wn_ref,'k^','MarkerSize',13,'MarkerFaceColor','k','Parent',Ax(1));
hold on
L2 = plot(rec00(:,1),imc00(:,1),'ko','MarkerFaceColor','k','MarkerSize',13, ...
    'LineWidth',.5,'Parent',Ax(1));
L4 = plot(real(S(:,1)'),imag(S(:,1)'),'pentagram','MarkerFaceColor','k', ...
    'MarkerEdgeColor','k','MarkerSize',13,'LineWidth',2);
set(gca,'FontSize',17);
xlim([-0.067 .067]);
ylim([.6 1.2]);
set(gca,'FontSize',21);
xlabel('\sigma_{1}/\omega_{r}')
ylabel('\omega_{1}/\omega_{r}')
grid on

Ax(2) = copyobj(Ax(1),gcf);
set(Ax(2),'Box', 'off','Visible','off','Position',get(Ax(1),'Position'))
delete(get(Ax(2),'Children'));
hold on

cc = jet(length(U));
plot(sigma_ref(:,1),wn_ref(:,1),'^','MarkerSize',13,'MarkerFaceColor', ...
    cc(20,:),'MarkerEdgeColor','k','HandleVisibility','off');
plot(sigma_ref(:,2),wn_ref(:,2),'^','MarkerSize',13,'MarkerFaceColor', ...
    cc(50,:),'MarkerEdgeColor','k','HandleVisibility','off');
plot(sigma_ref(:,3),wn_ref(:,3),'^','MarkerSize',13,'MarkerFaceColor', ...
    cc(80,:),'MarkerEdgeColor','k','HandleVisibility','off');

legend_entries{(length(U)/aux)} = ' ';
cont=1;
i_U = 1;
plot(rec00(:,i_U),imc00(:,i_U),'s','MarkerEdgeColor',cc(i_U,:),...
    'MarkerFaceColor',cc(i_U,:),'MarkerSize',10,'LineWidth',.5,'Parent',Ax(2));
hold on
legend_entries{cont} = sprintf('%0.2f',floor(Q(i_U)*10)/10);
cont=cont+1;
for i_U=aux:aux:length(U)
    plot(rec00(:,i_U),imc00(:,i_U),'s','MarkerEdgeColor',cc(i_U,:),...
        'MarkerFaceColor',cc(i_U,:),'MarkerSize',10,'LineWidth',.5,'Parent',Ax(2));
    legend_entries{cont} = sprintf('%0.2f',ceil(Q(i_U)*10)/10);
    cont=cont+1;
    hold on
end
leg = legend(Ax(2),legend_entries,'Orientation','vertical','location','EastOutside');

cont=1;
i_U = 1;
L3 = plot(real(S(:,i_U)'),imag(S(:,i_U)'),'pentagram','MarkerFaceColor',cc(i_U,:),...
    'MarkerSize',13,'MarkerEdgeColor','k','LineWidth',.5,'Parent',Ax(2),'HandleVisibility','off');
hold on
plot(rec00(:,i_U),imc00(:,i_U),'o','MarkerFaceColor',cc(i_U,:),'MarkerSize', ...
    13,'MarkerEdgeColor','k','LineWidth',.5,'Parent',Ax(1));
legend_entries{cont} = sprintf('%0.2f',floor(Q(i_U)*10)/10);
cont=cont+1;
for i_U=aux:aux:length(U)
    plot(real(S(:,i_U)'),imag(S(:,i_U)'),'pentagram','MarkerFaceColor',cc(i_U,:),...
        'MarkerSize',13,'MarkerEdgeColor','k','LineWidth',.5,'Parent',Ax(2), ...
        'HandleVisibility','off');
    hold on
    plot(rec00(:,i_U),imc00(:,i_U),'o','MarkerFaceColor',cc(i_U,:), ...
        'MarkerSize',13,'MarkerEdgeColor','k','LineWidth',.5,'Parent',Ax(1));
end

% legend(Ax(1),[L1 L2 L4],'Literature','SIMO','MIMO','location','northwest');
linkprop([Ax(1) Ax(2)], 'Position');
title(leg,'Q^{*}','fontweight','normal');
xlim([-0.067 .067]);
ylim([.6 1.2]);
xticks([-.06 -.03 0 .03 .06])
set(gca,'FontSize',21);
xlabel('\sigma_{1}/\omega_{r}')
ylabel('\omega_{1}/\omega_{r}')
grid on

figname = sprintf('Case%s_rLocus_1stMode',caso);
figures_save