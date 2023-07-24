%% Plot Impulse Responses for All Methods

% OKID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = -h./amp_pg;
u2 = aoa./amp_pt; 
y1 = cl;
y2 = cm;
u = [u1,u2]';
y = [y1,y2]';
addpath('functions');
[H,~] = OKID(y,u,r);

% plot OKID impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Clh = reshape(H(1,1,:),length(H),1);
Cla = reshape(H(1,2,:),length(H),1);
Cmh = reshape(H(2,1,:),length(H),1);
Cma = reshape(H(2,2,:),length(H),1);

% plot impulse responses in curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(linspace(0,0.1,length(H)+1),[0;Clh],'LineWidth',3)
ylabel({'$C_{\ell_h}$'}, 'interpreter', 'latex')
xlabel('Dimensionless Time')
ylim([-7e-4 13e-4])
yticks([-5e-4 0 5e-4 1e-3])
xticks([0 .02 .04 .06 .08 .1])
set(gca,'FontSize',21)
% figname = sprintf('OKIDImpResp11_r%d',r);
% dir = [fileparts(fileparts(pwd)),filesep,'3_output\r=' num2str(r) '\'];
% figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(linspace(0,0.1,length(H)+1),[0;Cla],'LineWidth',3)
ylabel({'$C_{\ell_\alpha}$'}, 'interpreter', 'latex')
xlabel('Dimensionless Time')
ylim([-3.5e-4 6.5e-4])
yticks([-2e-4 0 2e-4 4e-4 6e-4])
xticks([0 .02 .04 .06 .08 .1])
set(gca,'FontSize',21)
% figname = sprintf('OKIDImpResp12_r%d',r);
% figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(linspace(0,0.1,length(H)+1),[0;Cmh],'LineWidth',3)
ylabel({'$C_{m_h}$'}, 'interpreter', 'latex')
xlabel('Dimensionless Time')
ylim([-3.7e-4 2.2e-4])
yticks([-4e-4 -2e-4 0 2e-4])
xticks([0 .02 .04 .06 .08 .1])
set(gca,'FontSize',21)
% figname = sprintf('OKIDImpResp21_r%d',r);
% figures_save

figure
set(groot,'defaultAxesTickLabelInterpreter','latex');
plot(linspace(0,0.1,length(H)+1),[0;Cma],'LineWidth',3)
ylabel({'$C_{m_\alpha}$'}, 'interpreter', 'latex')
xlabel('Dimensionless Time')
ylim([-4e-4 2.5e-4])
yticks([-4e-4 -2e-4 0 2e-4])
xticks([0 .02 .04 .06 .08 .1])
set(gca,'FontSize',21)
% figname = sprintf('OKIDImpResp22_r%d',r);
% figures_save
