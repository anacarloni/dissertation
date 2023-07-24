clc
clearvars
close all

%% Figuras Resposta ao Impulso para Todos os ROMs

% Rank, r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=2:10

    % Respostas DS e US a partir do CFD, Baixa Amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idERA/yDS_lowAmp.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERA/yUS.mat'])

    % ROMs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idBPOD/sysBPOD_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERAOKID/sysERAOKID_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERA/sysERAlowAmp_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERA/sysERAmedAmp_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idERA/sysERAhighAmp_r' num2str(r) '.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_yVolt.mat'])

    % Calcular Respostas ao Impulso
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltat = 3e-3;
    tfinal = 300;
    t = 0:deltat:tfinal-deltat;
    tf = length(yDS_lowAmp)*deltat;
    y1 = permute(yDS_lowAmp,[3,1,2]);
    [y2,t2] = impulse(sysERAlowAmp,0:tf);
    [y3,t3] = impulse(sysBPOD,0:tf);
    [y4,t4] = impulse(sysERAOKID,0:tf);
    y5 = permute(yVolt,[3,1,2]);
    y2 = [zeros(1,2,2); y2];
    y3 = [zeros(1,2,2); y3];
    y4 = [zeros(1,2,2); y4];
    t2 = [t2; t2(end)+1];
    t3 = [t3; t3(end)+1];
    t4 = [t4; t4(end)+1];
    tmax = 0.1;
    tempo = t(t<tmax);
    t2n = t2.*deltat;
    t3n = t3.*deltat;
    t4n = t4.*deltat;

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg_lowAmp = str2double(amp_plunge);
    amp_pt_lowAmp = str2double(amp_pitch);
    amp_plunge = '500e-6';
    amp_pitch = '500e-4';
    amp_pg_medAmp = str2double(amp_plunge);
    amp_pt_medAmp = str2double(amp_pitch);
    amp_plunge = '1000e-6';
    amp_pitch = '1000e-4';
    amp_pg_highAmp = str2double(amp_plunge);
    amp_pt_highAmp = str2double(amp_pitch);

    % Resposta para Cada Amplitude usando Volterra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_ylow.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_ymed.mat'])
    load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_yhigh.mat'])

    % Figura - Resposta ao Impulso
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run('./figures/figures_impResp_st') % stairs
    run('./figures/figures_impResp')

    %% Plot da Resposta ao Impulso usando ERA para Diferentes Amplitudes

    % Calcular Respostas ao Impulso
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y1 = permute(yDS_lowAmp,[3,1,2]);
    [y2,t2] = impulse(sysERAlowAmp,0:length(yUS)*deltat);
    [y3,t3] = impulse(sysERAmedAmp,0:length(yUS)*deltat);
    [y4,t4] = impulse(sysERAhighAmp,0:length(yUS)*deltat);
    y2 = [zeros(1,2,2); y2];
    y3 = [zeros(1,2,2); y3];
    y4 = [zeros(1,2,2); y4];
    t2 = [t2; t2(end)+1];
    t3 = [t3; t3(end)+1];
    t4 = [t4; t4(end)+1];
    t2n = t2.*deltat;
    t3n = t3.*deltat;
    t4n = t4.*deltat;

    % Figura - Resposta ao Impulso
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run('./figures/figures_impResp_difAmp')
    close all
end