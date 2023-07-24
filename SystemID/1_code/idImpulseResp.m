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

    % Respostas DS e US a partir do CFD, Baixa Amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idERA/yUS.mat'])

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg = str2double(amp_plunge);
    amp_pt = str2double(amp_pitch);

    % Vetor de Tempo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = size(yUS,3);
    tempo = 0:N-1;
    deltat = 3e-3;
    tmax = 300;
    t = 0:deltat:tmax-deltat;

    % Resposta dos ROMs a um Sinal Arbitrário a partir do WF4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ht_imp = zeros(N,1);
    aoa_imp = zeros(N,1);
    ht_imp(2) = 1;
    aoa_imp(2) = 1;

    u = [ht_imp zeros(N,1)];
    yBPOD_pg = lsim(sysBPOD, u, tempo);
    yERAlowAmp_pg = lsim(sysERAlowAmp, u, tempo);
    yERAOKID_pg = lsim(sysERAOKID, u, tempo);

    u = [zeros(N,1) aoa_imp];
    yBPOD_pt = lsim(sysBPOD, u, tempo);
    yERAlowAmp_pt = lsim(sysERAlowAmp, u, tempo);
    yERAOKID_pt = lsim(sysERAOKID, u, tempo);

    % Figura - Resposta ao Impulso dos ROMs a partir do DS comparado com CFD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run('figures/figures_USresp')
    close all
end