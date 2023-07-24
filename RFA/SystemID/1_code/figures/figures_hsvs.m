%% Figuras dos Hankel Singular Values
mco = floor((length(yDS)/1000-1));
numInputs = 2;
numOutputs = 2;
YY = yDS(:,:,2:2*mco+2);
[~,~,~,~,hsvs2] = ERA(YY,mco,mco,numInputs,numOutputs,2);
[~,~,~,~,hsvs5] = ERA(YY,mco,mco,numInputs,numOutputs,5);
[~,~,~,~,hsvs6] = ERA(YY,mco,mco,numInputs,numOutputs,6);
[~,~,~,~,hsvs10] = ERA(YY,mco,mco,numInputs,numOutputs,10);

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
semilogy(hsvs2,'b*')
hold on, grid on
xlabel('$r$','Interpreter','latex')
ylabel('SV Magnitude')
ylim padded
yticks([1e-15 1e-10 1e-5])
set(gca,'FontSize',21)
figname = 'HankelSVs_magnitude';
dir = [fileparts(fileparts(pwd)),filesep,'3_output\hankelSV\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
semilogy(hsvs2,'b*')
hold on, grid on
xlabel('$r$','Interpreter','latex')
ylabel('SV Magnitude')
ylim([1e-6 1.5e-3])
xlim([0.4 17.9])
xticks([1 5 10 15 20])
yticks([1e-6 1e-5 1e-4 1e-3])
set(gca,'FontSize',21)
figname = 'zoom_HankelSVs_magnitude';
dir = [fileparts(fileparts(pwd)),filesep,'3_output\hankelSV\'];
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(0:length(hsvs2),[0; cumsum(hsvs2)/sum(hsvs2)],'b*')
% hold on
grid on
% plot(2,sum(hsvs2(1:2))/sum(hsvs2),'d','LineWidth',2,'MarkerSize',10)
% plot(5,sum(hsvs2(1:5))/sum(hsvs5),'s','LineWidth',2,'MarkerSize',10)
% plot(6,sum(hsvs2(1:6))/sum(hsvs6),'^','LineWidth',2,'MarkerSize',10)
% plot(10,sum(hsvs2(1:10))/sum(hsvs10),'o','LineWidth',2,'MarkerSize',10)
ylim padded
xlabel('$r$','Interpreter','latex')
ylabel('Energy')
set(gca,'FontSize',21)
% legend('', 'r=2', 'r=5', 'r=6', 'r=10','FontSize',21,'location','southeast')
figname = 'HankelSVs_energy';
figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(0:length(hsvs2),[0; cumsum(hsvs2)/sum(hsvs2)],'b*')
% hold on
grid on
% plot(2,sum(hsvs2(1:2))/sum(hsvs2),'d','LineWidth',2,'MarkerSize',10)
% plot(5,sum(hsvs2(1:5))/sum(hsvs5),'s','LineWidth',2,'MarkerSize',10)
% plot(6,sum(hsvs2(1:6))/sum(hsvs6),'^','LineWidth',2,'MarkerSize',10)
% plot(10,sum(hsvs2(1:10))/sum(hsvs10),'o','LineWidth',2,'MarkerSize',10)
ylim([.65 1.01])
xlim([0.4 17.9])
xticks([1 5 10 15 20])
xlabel('$r$','Interpreter','latex')
ylabel('Energy')
set(gca,'FontSize',21)
% legend('', 'r=2', 'r=5', 'r=6', 'r=10','FontSize',21,'location','southeast')
figname = 'zoom_HankelSVs_energy';
figures_save

% %% creating the zoom-in inset
% ax=axes;
% set(ax,'units','normalized','position',[0.3,0.3,0.5,0.5])
% box(ax,'on')
% plot(0:length(hsvs2),[0; cumsum(hsvs2)/sum(hsvs2)],'b*')
% set(ax,'xlim',[0 17.9],'ylim',[.65 1.01])
% xticks([1 5 10 15 20])
% grid on
