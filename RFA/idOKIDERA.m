clearvars
close all
clc

%% Identificar ROM usando ERA/OKID a partir da Resposta a Função de Walsh (WF)

% Rank, r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=2:10

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg = str2double(amp_plunge);
    amp_pt = str2double(amp_pitch);

    % Dados do Fort.14 (WF, Caso 13)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileID = fopen([fileparts(pwd),filesep,'0_data/Case13_Fort.14']);
    C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fileID);
    [t,ht,aoa,dcl,cd,dcm,cl,cm,cmex] = deal(C{:,1:9});
    N = size(t,1);
    Ts = t(2) - t(1);
    Fs = 1/Ts;

    % Variação dos Sinais de Entrada
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dht = -ht;
    daoa = aoa;
    for j = N:-1:2
        dht(j,1) = ht(j,1) - ht(j-1,1);
        daoa(j,1) = aoa(j,1) - aoa(j-1,1);
    end

    % Pressão Dinâmica
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho_inf = 1;
    alpha0 = 0;
    Mach = 0.8;
    U_inf = Mach*cos(alpha0);
    q = 0.5*rho_inf*U_inf^2;

    % Sinais de Entrada e Saída
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u1 = -ht./amp_pg;
    u2 = aoa./amp_pt;
    y1 = cl;
    y2 = cm;

    % Matriz de Resposta ao Impulso usando OKID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rERAOKID =  r;
    uRandom = [u1, u2]';
    yRandom = [y1, y2]';
    addpath('functions');
    [H,M] = OKID(yRandom,uRandom,rERAOKID);

    % Matrizes do ROM usando ERA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numInputs = 2;
    numOutputs = 2;
    mco = floor((length(H)-1)/2);
    [Ar,Br,Cr,Dr,HSVs] = ERA(H,mco,mco,numInputs,numOutputs,rERAOKID);
    sysERAOKID = ss(Ar,Br,Cr,Dr,-1);
    filename = sprintf('2_pipeline/idERAOKID/sysERAOKID_r%d.mat',r);
    save([fileparts(pwd),filesep,filename],'sysERAOKID');
end