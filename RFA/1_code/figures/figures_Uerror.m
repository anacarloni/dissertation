% Erro % entre U* para cada caso com n_beta polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(2:Npmax,array_er(:,1),'o-','MarkerSize',13,'LineWidth',2)
hold on
plot(2:Npmax,array_er(:,2),'+-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,3),'*-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,4),'d-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,5),'s-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,6),'^-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,7),'v-','MarkerSize',13,'LineWidth',2)
plot(2:Npmax,array_er(:,8),'<-','MarkerSize',13,'LineWidth',2)
xlim padded
ylim padded
xticks([1 2 3 4 5 6])
grid minor
xlabel('$n_\beta$','interpreter','latex')
ylabel('Error \%','interpreter','latex')
legend('Case 04','Case 05','Case 11','Case 12','Case 13','Case 16', ...
    'Case 17','Case 18','location','northoutside','Orientation','horizontal')
set(gca,'fontsize',21)
set(gcf, 'Position', get(0, 'Screensize'));
filename = '3_output/comparison/';
dir = [fileparts(pwd),filesep,filename];
figname = sprintf('Uerror_%s',form);
figures_save

%% creating the zoom-in inset
% ax=axes;
% set(ax,'units','normalized','position',[0.5,0.5,0.2,0.3])
% box(ax,'on')
% plot(1:Npmax,array_er(:,1),'o-','MarkerSize',13,'LineWidth',2,'parent',ax)
% hold on
% plot(1:Npmax,array_er(:,2),'+-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,3),'*-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,4),'d-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,5),'s-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,6),'^-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,7),'v-','MarkerSize',13,'LineWidth',2,'parent',ax)
% plot(1:Npmax,array_er(:,8),'<-','MarkerSize',13,'LineWidth',2,'parent',ax)
% set(ax,'xlim',[5.75 6.25],'ylim',[5 15])
% xticks([1 2 3 4 5 6])
% grid off