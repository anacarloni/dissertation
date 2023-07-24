% Plot norma L2 entre WF e RFA, para cada coeficiente aerodinâmico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(1:6,WF(1,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
% plot(2:6,RFAUnopt(1,2:end),'ks-','MarkerSize',13)
plot(1:6,RFA1st(1,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(1,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(1,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
grid off
ylim padded
xlim padded
xticks([1 2 3 4 5 6])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
set(gca,'FontSize',21)

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.58,0.3,0.3])
box(ax,'on')
plot(1:6,WF(1,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
plot(1:6,RFA1st(1,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(1,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(1,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
set(ax,'xlim',[3.8 6.2],'ylim',[.95 1.05])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
xticks([4 5 6])
set(gca,'FontSize',14) 
grid off

figname = sprintf('Case%s_L2norm_Clh',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(1:6,WF(2,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
% plot(2:6,RFAUnopt(2,2:end),'ks-','MarkerSize',13)
plot(1:6,RFA1st(2,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(2,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(2,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
grid off
if strcmp(caso,'16') || strcmp(caso,'17')
    ylim([min(WF(2,:))-0.2 max(RFA1st(:))+0.4])
else
    ylim padded
end
xlim padded
xticks([1 2 3 4 5 6])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
legend('Walsh function','First form','Second form','Third form')
set(gca,'FontSize',21)

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.34,0.3,0.3])
box(ax,'on')
plot(1:6,WF(2,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
plot(1:6,RFA1st(2,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(2,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(2,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
set(ax,'xlim',[3.8 6.2],'ylim',[.7 1.5])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
xticks([4 5 6])
set(gca,'FontSize',14)
grid off

figname = sprintf('Case%s_L2norm_Cla',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(1:6,WF(3,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
% plot(2:6,RFAUnopt(3,2:end),'ks-','MarkerSize',13)
plot(1:6,RFA1st(3,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(3,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(3,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
grid off
ylim padded
xlim padded
xticks([1 2 3 4 5 6])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
set(gca,'FontSize',21)

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.58,0.3,0.3])
box(ax,'on')
plot(1:6,WF(3,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
plot(1:6,RFA1st(3,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(3,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(3,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
set(ax,'xlim',[3.8 6.2],'ylim',[.4 .5])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
xticks([4 5 6])
set(gca,'FontSize',14)
grid off

figname = sprintf('Case%s_L2norm_Cmh',caso);
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(1:6,WF(4,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
% plot(2:6,RFAUnopt(4,2:end),'ks-','MarkerSize',13)
plot(1:6,RFA1st(4,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(4,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(4,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
grid off
ylim padded
xlim padded
xticks([1 2 3 4 5 6])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
set(gca,'FontSize',21)

% creating the zoom-in inset
ax=axes;
set(ax,'units','normalized','position',[0.55,0.59,0.3,0.3])
box(ax,'on')
plot(1:6,WF(4,:),'m*-','linewidth',2,'MarkerSize',13)
hold on
plot(1:6,RFA1st(4,:),'ro--','linewidth',2,'MarkerSize',13)
plot(1:6,RFA2nd(4,:),'b^-','linewidth',2,'MarkerSize',13)
plot(1,RFA2nd(4,1),'b^-','linewidth',2,'MarkerSize',13,'MarkerFaceColor','b')
set(ax,'xlim',[3.8 6.2],'ylim',[.47 .57])
xlabel('$n_\beta$','Interpreter','latex')
ylabel('Error')
xticks([4 5 6])
set(gca,'FontSize',14)
grid off

figname = sprintf('Case%s_L2norm_Cma',caso);
figures_save