% Comparação dos Resultados - Plano Complexo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = jet(length(k));
legend_entries{(floor(length(k)/aux))} = ' ';

f1 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
firstax = axes(f1); 
plot(real(kCLhfit),imag(kCLhfit),'k-','linewidth',2,'HandleVisibility','off');
hold all
plot(reclh,imclh,'k--','linewidth',2,'HandleVisibility','off')
for i=1:aux:length(k)
    plot(reclh(i),imclh(i),'x','Color',cc(i,:),...
        'MarkerSize',17,'LineWidth',2,'Parent',firstax);
    hold on
    legend_entries{(i-1)/aux+1} = sprintf('%0.2f',k(i));
end
leg = legend(legend_entries,'Orientation','vertical','location','northeast');
title(leg,'\kappa','FontSize',21)
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca, 'fontsize',21)
p1 = plot(reclh(id),imclh(id),'k^','MarkerSize',13,'HandleVisibility','off');
p2 = plot(real(CLhfit(id)),imag(CLhfit(id)),'k^','MarkerFaceColor','k', ...
    'MarkerSize',13,'HandleVisibility','off');
figname = sprintf('Case%s_IntPolFitComplexG11',caso);
figures_save

f1 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
firstax = axes(f1); 
p3 = plot(real(kCLafit),imag(kCLafit),'k-','linewidth',2,'HandleVisibility','off');
hold all
p4 = plot(recla,imcla,'k--','linewidth',2,'HandleVisibility','off');
for i=1:aux:length(k)
    plot(recla(i),imcla(i),'x','Color',cc(i,:),...
        'MarkerSize',17,'LineWidth',2, 'Parent', firstax);
    hold on
%     legend_entries{(i-1)/aux+1} = sprintf('%0.2f',k(i));
end
% leg = legend(firstax,legend_entries,'Orientation','vertical','location','bestoutside');
% title(leg,'\kappa','FontSize',21)
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca, 'fontsize',21)

secondax = copyobj(firstax, gcf);
linkprop([firstax secondax], 'Position');
p1 = plot(recla(id),imcla(id),'k^','MarkerSize',13, 'Parent', secondax);
p2 = plot(real(CLafit(id)),imag(CLafit(id)),'k^','MarkerFaceColor','k', ...
    'MarkerSize',13, 'Parent', secondax);
legend(secondax,[p3(1) p4(1)],'Data set','RFA fit','location','northeast');
figname = sprintf('Case%s_IntPolFitComplexG12',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(real(kCMhfit),imag(kCMhfit),'k-','linewidth',2,'HandleVisibility','off');
hold all
plot(recmh,imcmh,'k--','linewidth',2,'HandleVisibility','off')
for i=1:aux:length(k)
    plot(recmh(i),imcmh(i),'x','Color',cc(i,:),...
        'MarkerSize',17,'LineWidth',2);
    hold on
end
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca, 'fontsize',21)
p1 = plot(recmh(id),imcmh(id),'k^','MarkerSize',13);
p2 = plot(real(CMhfit(id)),imag(CMhfit(id)),'k^','MarkerFaceColor','k', ...
    'MarkerSize',13);
figname = sprintf('Case%s_IntPolFitComplexG21',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(real(kCMafit),imag(kCMafit),'k-','linewidth',2,'HandleVisibility','off');
hold all
plot(recma,imcma,'k--','linewidth',2,'HandleVisibility','off')
for i=1:aux:length(k)
    plot(recma(i),imcma(i),'x','Color',cc(i,:),...
        'MarkerSize',17,'LineWidth',2);
    hold on
end
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca, 'fontsize',21)
p1 = plot(recma(id),imcma(id),'k^','MarkerSize',13);
p2 = plot(real(CMafit(id)),imag(CMafit(id)),'k^','MarkerFaceColor','k', ...
    'MarkerSize',13);
figname = sprintf('Case%s_IntPolFitComplexG22',caso);
figures_save
