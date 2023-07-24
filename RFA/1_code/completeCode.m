clc
clearvars
close all

%% Funções de Transferência no Domínio da Frequência
casos = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09'; '10'; '11';
    '12'; '13'; '14'; '15'; '16'; '17'; '18'];
Npmax = 6;
for i=1:length(casos)
    caso = casos(i,:);
    switch caso
        case {'01','02'}
            inputsignal = 'WF1';
            nblocks = 2;
        case {'03','04','05'}
            inputsignal = 'WF2';
            nblocks = 3;
        case {'06','07','08'}
            inputsignal = 'WF3';
            nblocks = 3;
        case {'09','10','11','12','13'}
            inputsignal = 'WF4';
            nblocks = 4;
        case {'14','15','16','17','18'}
            inputsignal = 'WF5';
            nblocks = 4;
    end
    output_spectral
end

%% Primeira Forma dos RFAs - Polos Não Otimizados
casos = ['04'; '05'; '11'; '12'; '13'; '16'; '17'; '18'];
Npmax = 6;

% Polos Não Otimizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        unoptimized
    end
end

% Análise de Estabilidade Aeroelástica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        stability_analysis_unoptimized
    end
end

%% Primeira Forma dos RFAs - Polos Otimizados

% Polos Otimizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        optimized_1stForm
    end
end

% Análise de Estabilidade Aeroelástica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        stability_analysis_1stForm
    end
end

%% Segunda Forma dos RFAs - Polos Otimizados

% Polos Otimizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        optimized_2ndForm
    end
end

% Análise de Estabilidade Aeroelástica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(casos)
    caso = casos(i,:);
    for Np=1:Npmax
        stability_analysis_2ndForm
    end
end

%% Tabelas de Resultados
output_table