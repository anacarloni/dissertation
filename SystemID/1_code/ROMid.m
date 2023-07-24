clearvars
close all
clc
addpath('figures')

%% Problema Aeroelástico

% Parâmetros do Sistema Aeroelástico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = .5;         % semi-chord length
xam = 1.8;      % distance from the EA to the CM, normalized by semichord
xea = 1/2;      % EA position used to calculate aerodynamic coefficients, semi-chord
ah = -2;        % distance from midchord to the EA, normalized by semichord
ra = 1.865;     % airfoil dimensionless radius of gyration about the EA
wa = 100;       % uncoupled natural circular frequency of the pitch DOF, rad/s
wh = 100;       % uncoupled natural circular frequency of the plunge DOF, rad/s
wr = wa;        % reference circular frequency, rad/s
mu = 60;        % mass ratio, mu = m/(pi*rho_inf*b^2)
Mach = .8;      % freestream Mach number
ndof = 2;       % number of DOFs
rho_inf = 1;    % freestream density
m = mu*pi*rho_inf*b^2; % airfoil mass, normalized by reference length

% Pressão Dinâmica Característica, Q*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qmin = 0.1;
Qmax = 1;
NQ = 10;
Q = linspace(Qmin,Qmax,NQ);
Qref = 0.5;
Uref = sqrt(Qref*mu);

%% Excitação Simultânea usando Funções de Walsh (WF)

% Resposta do CFD a WF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casos = ['02'; '05'; '08'; '13'; '18']; % cada conjunto de WF
for v=1:length(casos)
    caso = casos(v,:);
    file = sprintf('0_data/Case%s_Fort.14',caso);
    fileID = fopen([fileparts(pwd),filesep,file]);
    C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fileID);
    [t,h,alpha,dcl,~,dcm,cl,cm,~] = deal(C{:,1:9});
    N = size(t,1);
    deltat = t(2) - t(1);

    % Forças Aerodinâmicas Generalizadas (GAFs) normalizadas pela pressão dinâmica
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fa = [(-2/m).*cl (4/m).*cm]./wr^2;

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg = str2double(amp_plunge);
    amp_pt = str2double(amp_pitch);


    %% Identificação do ROM Aerodinâmico

    % OKID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run('./figures/figures_OKIDimpResp'); % reproduce impulse responses
    u1 = -h;
    u2 = -alpha;
    y1 = Fa(:,1);
    y2 = Fa(:,2);
    u = [u1,u2]';
    y = [y1,y2]';
    addpath('functions');

    % Método de Identificação do ROM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods = ["OKIDERA"; "BPOD"];
    rank = 2:10;
    Qflutter = zeros(length(rank),1);
    Uflutter = zeros(length(rank),1);
    error = zeros(length(rank),1);
    for l=1:length(methods)
        for ind=1:length(rank)
            r = rank(ind);
            [H,~] = OKID(y,u,r);

            % Coeficientes Aerodinâmicos da Matriz H
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Clh = reshape(H(1,1,:),length(H),1);
            Cla = reshape(H(1,2,:),length(H),1);
            Cmh = reshape(H(2,1,:),length(H),1);
            Cma = reshape(H(2,2,:),length(H),1);

            % Corrigindo os Coeficientes (Posição do Eixo Elástico)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            arm = -ah-xea;
            Cma = Cma - arm*(Cla + Cmh - arm*Clh);
            Cla = Cla - arm*Clh;
            Cmh = Cmh - arm*Clh;
            for i=1:length(H)
                H(:,:,i) = [Clh(i) Cla(i); Cmh(i) Cma(i)];
            end

            % Identificação do ROM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idmethod = methods(l,:);
            switch idmethod
                case {'OKIDERA','BPOD'}
                    numInputs = 2;
                    numOutputs = 2;
                    mco = floor((length(H)-1)/2);
                    [A,B,C,D,HSVs] = ERA(H,mco,mco,numInputs,numOutputs,r);
                    aerosys = ss(A,B,C,D,-1);
                    filename = sprintf('2_pipeline/ROMid/Case%s/sysERA_r%d.mat', ...
                        caso,r);
                    save([fileparts(pwd),filesep,filename],'aerosys');

                    if(strcmp(idmethod,'BPOD'))
                        sysERA = aerosys;
                        [~,~,xFull] = impulse(sysERA,0:1:(r*5)+1);
                        sysAdjoint = ss(sysERA.A',sysERA.C',sysERA.B', ...
                            sysERA.D',-1);
                        [yAdjoint,~,xAdjoint] = impulse(sysAdjoint, ...
                            0:1:(r*5)+1);
                        HankelOC = [];
                        for i=2:size(xAdjoint,1)
                            Hrow = [];
                            for j=2:size(xFull,1)
                                Ystar = permute(squeeze(xAdjoint(i,:,:)), ...
                                    [2 1]);
                                MarkovParameter = Ystar*squeeze(xFull(j,:,:));
                                Hrow = [Hrow MarkovParameter];
                            end
                            HankelOC = [HankelOC; Hrow];
                        end
                        [U,Sig,V] = svd(HankelOC);
                        Xdata = [];
                        Ydata = [];
                        for i=2:size(xFull,1)
                            Xdata = [Xdata squeeze(xFull(i,:,:))];
                            Ydata = [Ydata squeeze(xAdjoint(i,:,:))];
                        end
                        Phi = Xdata*V*Sig^(-1/2);
                        Psi = Ydata*U*Sig^(-1/2);
                        A = Psi(:,1:r)'*sysERA.a*Phi(:,1:r);
                        B = Psi(:,1:r)'*sysERA.b;
                        C = sysERA.c*Phi(:,1:r);
                        D = sysERA.d;
                        aerosys = ss(A,B,C,D,-1);
                        filename = sprintf('2_pipeline/ROMid/Case%s/sysBPOD_r%d.mat', ...
                            caso,r);
                        save([fileparts(pwd),filesep,filename],'aerosys');
                    end

            end
            Aa = aerosys.A;
            Ba = aerosys.B;
            Ca = aerosys.C;
            Da = aerosys.D;

            %% ROM Estrutural
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M = [1 xam; xam ra^2];
            K = [((wh/wr)^2) 0 ; 0 (ra*wa/wr)^2];
            As = [zeros(ndof) eye(ndof); -M\K zeros(ndof)];
            Bs = [zeros(ndof); inv(M)];
            Cs = [eye(ndof) zeros(ndof)];
            strucsys = ss(As,Bs,Cs,[]);

            %% Análise de Estabilidade Aeroelástica

            % Alocação de Memória
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P = zeros(2*ndof+r,NQ);
            Z = zeros(2*ndof+r,NQ);
            stab = zeros(1,NQ);

            q = .5*mu*rho_inf*Q*(b*wr)^2;
            figure;
            set(groot,'defaultAxesTickLabelInterpreter','latex');
            for i=1:NQ
                % Modelo Estrutural no Tempo Discreto
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Uast(i) = sqrt(mu*Q(i));
                deltatbar = 2*Mach*deltat/Uast(i);
                structdsys = c2d(strucsys,deltatbar);
                barAs = structdsys.A;
                barBs = structdsys.B;

                % Cálculo de barAs e barBs de forma explícita
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % fun = @(x) expm(As.*x);
                % barAs = expm(As.*deltat);
                % barBs = integral(fun,0,deltat,'ArrayValued',1)*Bs;

                % Modelo Aeroelástico no Tempo Discreto
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                A = [barAs+q(i).*barBs*Da*Cs q(i).*barBs*Ca; Ba*Cs Aa];
                C = [Cs zeros(ndof,r)];
                sys_ae = ss(A,zeros(length(A),1),C,zeros(2,1),deltat);

                % Análise de Estabilidade
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % [p,Z] = pzmap(sys_ae);
                stab(i) = isstable(sys_ae);

                % Lugar das Raízes no Plano Z Complexo
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hz = pzplot(sys_ae);
                drawnow
                hold all
                hm = findobj(gca, 'Type', 'Line');
                hm(3).MarkerSize = 15;
                hm(3).LineWidth=2;
                legend_entries{(i)} = sprintf('%0.2f',Q(i));
            end
            p = getoptions(hz);
            p.Title.String = '';
            p.XLabel.FontSize = 21;
            p.YLabel.FontSize = 21;
            p.TickLabel.FontSize = 10;
            p.TickLabel.Color = [0 0 0];
            setoptions(hz,p);
            grid on
            set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 10], ...
                'PaperUnits', 'Inches', 'PaperSize', [8, 10]);
            leg = legend(legend_entries,'Orientation','vertical','location', ...
                'northeast');
            title(leg,'Q^{*}','fontweight','normal');

            % Círculo Unitário
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             th = 0:pi/50000:2*pi;
            %             rad = 1;
            %             xunit = rad * cos(th);
            %             yunit = rad * sin(th);
            %             plot(xunit, yunit,'k:','LineWidth',3,'HandleVisibility','off');
            %             xlim([-1.05 1.05])
            %             hz.getaxes.FontSize = 21;

            % Figura - Lugar das Raízes
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figname = sprintf('RootLoc_Case%s_%s_r%s',caso,idmethod, ...
                num2str(r));
            dir = [fileparts(pwd),filesep,'3_output\rootLocus\'];
            figures_save

            % Ponto de Início do Flutter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            array = Q(stab==0);
            disp('Flutter onset point (Q*):')
            if isempty(array)
                Qflutter(ind) = 0;
                Uflutter(ind) = 0;
                error(ind) = 0;
                disp('Stable');
            else
                Qflutter(ind) = array(1);
                Uflutter(ind) = sqrt(Qflutter(ind)*mu);
                error(ind) = (abs(Uflutter(ind)-Uref)./Uref)*100;
                disp(Qflutter(ind));
            end
            close all
        end

        format bank
        Uast = sqrt(mu*Q);
        legend(num2str(Q(1)),num2str(Q(2)),num2str(Q(3)),num2str(Q(4)), ...
            num2str(Q(5)),num2str(Q(6)),num2str(Q(7)), ...
            num2str(Q(8)),num2str(Q(9)),num2str(Q(10)));

        % Arquivos de Saída
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out = [rank; Qflutter; Uflutter; error];
        fout = fopen(sprintf('Case%s_%s_flutter.dat',caso,idmethod),'w');
        fprintf(fout,'VARIABLES = "Rank" "Q" "U" "Error"\n');
        fprintf(fout,'NQ=%4d\n',NQ);
        fprintf(fout,'%d %11.7E %11.7E %11.7E\n',out);
        fclose(fout);
        filename = [fileparts(pwd),filesep,'3_output/rootLocus/'];
        movefile('*.dat',filename)
        clc
    end
end