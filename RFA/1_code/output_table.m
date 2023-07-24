% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
close all
format shortG
addpath('figures');

%% Caso Específico para Todos os Números de Polos

% Casos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casos = ['04'; '05'; '11'; '12'; '13'; '16'; '17'; '18'];

% Quantidade de Polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npmin = 1;
Npmax = 6;
Np = (1:Npmax)';

% Formulação do RFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option = ["unopt" "1st" "2nd"]';
for j=1:length(option)
    if strcmp(option(j),'unopt')
        form = 'unoptimized';
    elseif strcmp(option(j),'1st')
        form = 'optimized_1stForm';
    else
        form = 'optimized_2ndForm';
    end

    % Velocidade de Flutter de Referência (Rausch, Batina and Yang, 1990)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu = 60;
    Qref = 0.5;
    Uref = sqrt(Qref*mu);

    % Alocação de Memória
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = zeros(Npmax,1);
    U = zeros(Npmax,1);
    s1 = zeros(Npmax,1);
    w1 = zeros(Npmax,1);
    s2 = zeros(Npmax,1);
    w2 = zeros(Npmax,1);
    error = zeros(Npmax,1);
    array_er = zeros(Npmax-1,size(casos,1));

    % Leitura dos Dados
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for a = 1:size(casos,1)
        caso = casos(a,:);
        Np = (1:Npmax)';
        for i = 1:Npmax
            filename = sprintf('2_pipeline/stability_analysis/%s/%d_poles/Case%s_flutter.dat',...
                form,Np(i),caso);
            if isfile([fileparts(pwd),filesep,filename])
                fileID = fopen([fileparts(pwd),filesep,filename]);
                C = textscan(fileID,'%f %f %f %f %f %f','HeaderLines',3);
                fclose(fileID);
                [Q(i),U(i),s1(i),w1(i),s2(i),w2(i)] = deal(C{:,:});
                error(i) = (abs(U(i)-Uref)./Uref)*100;
            end
        end

        % Tabela com os Resultados
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Np=Np(2:Npmax); s1=s1(2:Npmax); s2=s2(2:Npmax); w1=w1(2:Npmax);
        w2=w2(2:Npmax); Q=Q(2:Npmax); U=U(2:Npmax); error=error(2:Npmax);
        array_er(:,a) = error;
        tabela = table(Np,s1,w1,s2,w2,Q,U,error);
        disp(tabela);
        dir = [fileparts(pwd),filesep,sprintf('3_output/%s/',form)];
        writetable(tabela,sprintf('Case%s_table.dat',caso),'WriteRowNames', ...
            true, ...
            'Delimiter',' ');
        status = movefile(sprintf('Case%s_table.dat',caso),dir,'f');
    end
    figures_Uerror
end

%% Número de Polos Específico para Todos os Casos

% Casos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casos = ['04'; '05'; '11'; '12'; '13'; '16'; '17'; '18'];

% Formulação do RFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option = ["unopt" "1st" "2nd"]';
for j=1:length(option)
    if strcmp(option(j),'unopt')
        form = 'unoptimized';
    elseif strcmp(option(j),'1st')
        form = 'optimized_1stForm';
    else
        form = 'optimized_2ndForm';
    end

    % Alocação de Memória
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dim = size(casos,1);
    Q = zeros(dim,1);
    U = zeros(dim,1);
    s1 = zeros(dim,1);
    w1 = zeros(dim,1);
    s2 = zeros(dim,1);
    w2 = zeros(dim,1);
    error = zeros(dim,1);

    % Leitura dos Dados
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Np = (1:Npmax)';
    for a = 1:Npmax
        for i = 1:size(casos,1)
            caso = casos(i,:);
            filename = sprintf('2_pipeline/stability_analysis/%s/%d_poles/Case%s_flutter.dat',...
                form,Np(a),caso);
            if isfile([fileparts(pwd),filesep,filename])
                fileID = fopen([fileparts(pwd),filesep,filename]);
                C = textscan(fileID,'%f %f %f %f %f %f','HeaderLines',3);
                fclose(fileID);
                [Q(i),U(i),s1(i),w1(i),s2(i),w2(i)] = deal(C{:,:});
                error(i) = (abs(U(i)-Uref)./Uref)*100;
            else
                continue;
            end
        end

        % Tabela com os Resultados
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Case = casos;
        tabela = table(Case,s1,w1,s2,w2,Q,U,error);
        disp(tabela);
        dir = [fileparts(pwd),filesep,sprintf('3_output/%s/%s_poles/',form,...
            num2str(Np(a)))];
        writetable(tabela,'table.dat','WriteRowNames',true,'Delimiter',' ');
        status = movefile('table.dat',dir,'f');
    end
end

%% Norma L2 Comparado com Caso 00 (SIMO)

% Casos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casos = ['04'; '05'; '11'; '12'; '13'; '16'; '17'; '18'];

% Número de polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npmax = 6;

for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        % Leitura dos Dados de Norma L2 (Unoptimized)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        form = 'unoptimized';
        filename = sprintf('2_pipeline/%s/%d_poles/Case%s_L2norm.dat',...
            form,Np,caso);
        if isfile([fileparts(pwd),filesep,filename])
            fileID = fopen([fileparts(pwd),filesep,filename]);
            C = textscan(fileID,'%f %f %f %f','HeaderLines',1);
            fclose(fileID);
            [ClhUnopt,ClaUnopt,CmhUnopt,CmaUnopt] = deal(C{:,:});

            % Tabela com os Resultados
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tabela = array2table([ClhUnopt,ClaUnopt,CmhUnopt,CmaUnopt], ...
                'RowNames', {'WF', 'RFA'});
            tabela.Properties.VariableNames = ["Clh","Cla","Cmh","Cma"];
            % disp(tabela);
            dir = [fileparts(pwd),filesep,sprintf('3_output/%s/%s_poles/', ...
                form,num2str(Np))];
            writetable(tabela,sprintf('Case%s_tableL2norm.dat',caso), ...
                'WriteRowNames', ...
                true,'Delimiter',' ');
            status = movefile(sprintf('Case%s_tableL2norm.dat',caso),dir, ...
                'f');
        end

        % Leitura dos Dados de Norma L2 (1st Form)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        form = 'optimized_1stForm';
        filename = sprintf('2_pipeline/%s/%d_poles/Case%s_L2norm.dat',...
            form,Np,caso);
        if isfile([fileparts(pwd),filesep,filename])
            fileID = fopen([fileparts(pwd),filesep,filename]);
            C = textscan(fileID,'%f %f %f %f','HeaderLines',1);
            fclose(fileID);
            [Clh1st,Cla1st,Cmh1st,Cma1st] = deal(C{:,:});

            % Tabela com os Resultados
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tabela = array2table([Clh1st,Cla1st,Cmh1st,Cma1st],'RowNames', ...
                {'WF', 'RFA'});
            tabela.Properties.VariableNames = ["Clh","Cla","Cmh","Cma"];
            % disp(tabela);
            dir = [fileparts(pwd),filesep,sprintf('3_output/%s/%s_poles/', ...
                form,num2str(Np))];
            writetable(tabela,sprintf('Case%s_tableL2norm.dat',caso), ...
                'WriteRowNames', ...
                true,'Delimiter',' ');
            status = movefile(sprintf('Case%s_tableL2norm.dat',caso),dir, ...
                'f');
        end

        % Leitura dos Dados de Norma L2 (2nd Form)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        form = 'optimized_2ndForm';
        filename = sprintf('2_pipeline/%s/%d_poles/Case%s_L2norm.dat',...
            form,Np,caso);
        if isfile([fileparts(pwd),filesep,filename])
            fileID = fopen([fileparts(pwd),filesep,filename]);
            C = textscan(fileID,'%f %f %f %f','HeaderLines',1);
            fclose(fileID);
            [Clh2nd,Cla2nd,Cmh2nd,Cma2nd] = deal(C{:,:});

            % Tabela com os Resultados
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tabela = array2table([Clh2nd,Cla2nd,Cmh2nd,Cma2nd], 'RowNames', ...
                {'WF', 'RFA'});
            tabela.Properties.VariableNames = ["Clh","Cla","Cmh","Cma"];
            % disp(tabela);
            dir = [fileparts(pwd),filesep,sprintf('3_output/%s/%s_poles/', ...
                form,num2str(Np))];
            writetable(tabela,sprintf('Case%s_tableL2norm.dat',caso), ...
                'WriteRowNames', ...
                true,'Delimiter',' ');
            status = movefile(sprintf('Case%s_tableL2norm.dat',caso),dir, ...
                'f');
        end

        % Plot norma L2 entre WF e RFA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         WF = [Clh1st(1),Cla1st(1),Cmh1st(1),Cma1st(1)];
        %         RFA1st = [Clh1st(2),Cla1st(2),Cmh1st(2),Cma1st(2)];
        %         RFA2nd = [Clh2nd(2),Cla2nd(2),Cmh2nd(2),Cma2nd(2)];
        %         addpath('figures');
        %         dir = [fileparts(pwd),filesep,'3_output/comparison/' num2str(Np) '_poles\'];
        %         figures_L2norm_WFandRFA
        %         close all

        % Plot norma L2 entre WF e RFA para cada coeficiente aerodinâmico
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WF(:,Np) = [Clh1st(1);Cla1st(1);Cmh1st(1);Cma1st(1)];
        RFAUnopt(:,Np) = [ClhUnopt(2);ClaUnopt(2);CmhUnopt(2);CmaUnopt(2)];
        RFA1st(:,Np) = [Clh1st(2);Cla1st(2);Cmh1st(2);Cma1st(2)];
        RFA2nd(:,Np) = [Clh2nd(2);Cla2nd(2);Cmh2nd(2);Cma2nd(2)];
    end
    addpath('figures');
    dir = [fileparts(pwd),filesep,'3_output/comparison/'];
    figures_L2norm_coefs
    close all
end