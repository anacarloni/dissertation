% Block/sample/overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('functions')

figure
plot(linspace(1,N,1001),[0 h_inp(1:100:end)'./(abs(max(h_inp)))],'b', ...
    'LineWidth',2,'HandleVisibility','off')
hold on
plot(0*ones(1,10), linspace(-2,2.1,10),'r--','LineWidth',2)
plot((N/3)*ones(1,10), linspace(-2,2.1,10),'r--','LineWidth',2)
plot((2*N/3)*ones(1,10), linspace(-2,2.1,10),'r--','LineWidth',2)
plot(N*ones(1,10), linspace(-2,2.1,10),'r--','LineWidth',2)
ylabel('h/|h_{max}|')
xlabel('Time Step')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',21);

% figname = sprintf('input_%s_Case04_blocks',inputsignal);
% figures_save

figure
plot(linspace(1,N,1001),[0 h_inp(1:100:end)'./(abs(max(h_inp)))],'b', ...
    'LineWidth',2,'HandleVisibility','off')
hold on
plot(0*ones(1,10), linspace(-2,2,10),'g--','LineWidth',2)
plot((2*N/3)*ones(1,10), linspace(-2,2,10),'g--','LineWidth',2,'HandleVisibility','off')
plot((N/3)*ones(1,10), linspace(-2,2,10),'m-.','LineWidth',2)
plot(N*ones(1,10), linspace(-2,2,10),'m-.','LineWidth',2,'HandleVisibility','off')
ylabel('h/|h_{max}|')
xlabel('Time Step')
legend('Sample 1','Sample 2')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',21);

% figname = sprintf('input_%s_Case04_samples',inputsignal);
% figures_save

% Figuras para gerar animação no TeX - Case 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
cyan = [0.3010 0.7450 0.9330];
red = [0.6350 0.0780 0.1840];

figure
plot(linspace(1,N,1001),[0 h_inp(1:100:end)']./(abs(max(h_inp))),'b','LineWidth',2)
ylabel('h/|h_{max}|')
xlabel('Time Step')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',17);

% figname = sprintf('input_%s_Case11',inputsignal);
% figures_save

figure
plot(linspace(1,N,1001),[0 h_inp(1:100:end)'./(abs(max(h_inp)))],'b','LineWidth',2)
ylabel('h/|h_{max}|')
xlabel('Time Step')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',17);
ShadePlot(0, N/4, yellow)
ShadePlot(N/4+1, 2*N/4, purple)
ShadePlot(2*N/4+1, 3*N/4, green)
ShadePlot(3*N/4+1, 4*N/4, cyan)

% figname = sprintf('input_%s_Case11_blocks',inputsignal);
% figures_save

figure
plot(linspace(1,N,1000),h_inp(1:100:end)'./(abs(max(h_inp))),'b','LineWidth',2)
ylabel('h/|h_{max}|')
xlabel('Time Step')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',17);
% ShadePlot(0, 2*N/4, yellow)
% ShadePlot(N/4+1, 3*N/4, purple)
ShadePlot(2*N/4+1, 4*N/4, green)

% figname = sprintf('input_%s_Case11_samples3',inputsignal);
% figures_save

figure
plot(linspace(1,N,1000),h_inp(1:100:end)'./(abs(max(h_inp))),'b','LineWidth',2)
ylabel('h/|h_{max}|')
xlabel('Time Step')
ylim([-1.5 1.5])
xlim([0 N])
grid on
set(gca,'FontSize',17);
% ShadePlot(0, N/4, yellow)
ShadePlot(N/4+1, 2*N/4, purple)
ShadePlot(3*N/4+1, 4*N/4, green)
ShadePlot(2*N/4+1, 3*N/4, 'r') % Overlap

% figname = sprintf('input_%s_Case11_overlap2',inputsignal);
% figures_save