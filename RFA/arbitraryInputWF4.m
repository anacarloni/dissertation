clearvars
close all
clc

%% Modelos de Ordem Reduzida (ROMs) - ERA, BPOD, OKID/ERA

% Rank, r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=2:10

    % ROMs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idBPOD/sysBPOD_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERAOKID/sysERAOKID_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERA/sysERAlowAmp_r' num2str(r) '.mat'])

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg_lowAmp = str2double(amp_plunge);
    amp_pt_lowAmp = str2double(amp_pitch);

    % Dados do Case13_Fort.14 (WF4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = '0_data/Case13_Fort.14';
    fileID = fopen([fileparts(pwd),filesep,filename]);
    C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fileID);
    [t,ht,aoa,dcl,~,dcm,cl,cm,~] = deal(C{:,1:9});
    N = size(t,1);
    tempo = 0:N-1;

    % Variação dos Sinais de Entrada
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dht = -ht;
    daoa = aoa;
    for j = 2:N
        dht(j,1) = ht(j,1) - ht(j-1,1);
        daoa(j,1) = aoa(j,1) - aoa(j-1,1);
    end
    dht = flip(dht,1);
    daoa = flip(daoa,1);
    dht = [zeros(3,1); dht(1:end-3)];
    daoa = [zeros(3,1); daoa(1:end-3)];

    % Convenção de Sinal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcl = -dcl;
    dcm = -dcm;
    dcl(1) = 0;
    dcm(1) = 0;

    % Resposta a Função de Walsh WF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yWF4 = [dcl dcm];

    % Resposta dos ROMs ao WF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = [dht./amp_pg_lowAmp daoa./amp_pt_lowAmp];
    yBPOD = lsim(sysBPOD,u,tempo);
    yERAlowAmp = lsim(sysERAlowAmp,u,tempo);
    yERAOKID = lsim(sysERAOKID,u,tempo);

    % Figura - Resposta ao Impulso dos ROMs a partir do DS comparado com CFD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run('figures/figures_WF4resp')
    close all
end