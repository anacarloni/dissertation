% Comparação dos Resultados - Plano Complexo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = jet(length(k));
legend_entries{(floor(length(k)/aux))} = ' ';

f1 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
firstax = axes(f1); 
plot(real(kCLhfit),imag(kCLhfit),'k-','HandleVisibility','off'); 
hold on
plot(reclh,imclh,'k--','HandleVisibility','off')
for i=1:aux:length(k)
    plot(reclh(i),imclh(i),'x','Color',cc(i,:),...
        'MarkerSize',13,'LineWidth',2,'Parent',firstax);
    hold on
    legend_entries{(i-1)/aux+1} = sprintf('%0.2f',k(i));
end
leg = legend(legend_entries,'Orientation','vertical','location','northeast');
title(leg,'\kappa','FontSize',21)
ylabel('Imag')
xlabel('Real')
xlim padded
ylim padded
grid on
set(gca,'fontsize',21)
figname = sprintf('Case%s_IntPolComplexG11',caso);
figures_save

f1 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
firstax = axes(f1); 
p1 = plot(real(kCLafit),imag(kCLafit),'k-','HandleVisibility','off');
hold on
p2 = plot(recla,imcla,'k--','HandleVisibility','off');
for i=1:aux:length(k)
    plot(recla(i),imcla(i),'x','Color',cc(i,:),...
        'MarkerSize',13,'LineWidth',2, 'Parent', firstax);
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
set(gca,'fontsize',21)

secondax = copyobj(firstax, gcf);
linkprop([firstax secondax], 'Position');
legend(secondax,[p1 p2],'Data set','RFA fit','location','northeast');

figname = sprintf('Case%s_IntPolComplexG12',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(real(kCMhfit),imag(kCMhfit),'k-','HandleVisibility','off');
hold on
plot(recmh,imcmh,'k--','HandleVisibility','off')
for i=1:aux:length(k)
    plot(recmh(i),imcmh(i),'x','Color',cc(i,:),...
        'MarkerSize',13,'LineWidth',2);
    hold on
end
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca,'fontsize',21)
figname = sprintf('Case%s_IntPolComplexG21',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(real(kCMafit),imag(kCMafit),'k-','HandleVisibility','off');
hold on
plot(recma,imcma,'k--','HandleVisibility','off')
for i=1:aux:length(k)
    plot(recma(i),imcma(i),'x','Color',cc(i,:),...
        'MarkerSize',13,'LineWidth',2);
    hold on
end
xlim padded
ylim padded
ylabel('Imag')
xlabel('Real')
grid on
set(gca,'fontsize',21)
figname = sprintf('Case%s_IntPolComplexG22',caso);
figures_save

% lg  = legend(legend_entries,'Orientation','vertical'); 
% set(gcf, 'WindowState', 'maximized');
% lg.Position = [.96 0.54 0 0];
% lg.Title.String = '\kappa';
