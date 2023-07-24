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

    % Pole-Zero Map dos ROMs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure; pzmap(sysBPOD);
    % figure; pzmap(sysERAOKID);
    % figure; pzmap(sysERAlowAmp);

    %% Resposta do CFD ao Discrete Step (DS)

    % Amplitudes Máximas em Plunge e Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amp_plunge = '1e-6';
    amp_pitch = '1e-4';
    amp_pg_lowAmp = str2double(amp_plunge);
    amp_pt_lowAmp = str2double(amp_pitch);

    % Dados do Fort.14 (DS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename_pg = sprintf('0_data/DSplunge_%sc_Fort.14',amp_plunge);
    filename_pt = sprintf('0_data/DSpitch_%s_Fort.14',amp_pitch);
    fileID1 = fopen([fileparts(pwd),filesep,filename_pg]);
    fileID2 = fopen([fileparts(pwd),filesep,filename_pt]);
    C1 = textscan(fileID1,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
    C2 = textscan(fileID2,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
    fclose(fileID1);
    fclose(fileID2);
    [~,ht,~,dcl_pg,~,dcm_pg,cl_pg,cm_pg,~] = deal(C1{:,1:9});
    [t,~,aoa,dcl_pt,~,dcm_pt,cl_pt,cm_pt,~] = deal(C2{:,1:9});
    N = size(t,1);
    tempo = 0:N-1;

    % Convenção de Sinal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcl_pg = -dcl_pg;
    dcm_pg = -dcm_pg;

    % Variação dos Sinais de Entrada
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dht = -ht;
    daoa = aoa;
    for j = N:-1:2
        dht(j,1) = ht(j,1)-ht(j-1,1);
        daoa(j,1) = aoa(j,1)-aoa(j-1,1);
    end
    dht = [dht(2:end); dht(end)];
    daoa = [daoa(2:end); daoa(end)];

    % Resposta ao Impulso (DS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yDS = zeros(2,2,N);
    for i=2:N
        yDS(:,:,i-1) = [dcl_pg(i) dcl_pt(i); dcm_pg(i) dcm_pt(i)];
    end
    yDS(:,:,N) = yDS(:,:,end);

    %% Resposta dos ROMs response a Unit Sample (US) a partir do DS

    % Resposta dos ROMs para US Prescrito em Plunge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = [dht./amp_pg_lowAmp, zeros(N,1)];
    yBPOD_pg = lsim(sysBPOD,u,tempo);
    yERAlowAmp_pg = lsim(sysERAlowAmp,u,tempo);
    yERAOKID_pg = lsim(sysERAOKID,u,tempo);

    % Resposta dos ROMs para US Prescrito em Pitch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = [zeros(N,1), daoa./amp_pt_lowAmp];
    yBPOD_pt = lsim(sysBPOD,u,tempo);
    yERAlowAmp_pt = lsim(sysERAlowAmp,u,tempo);
    yERAOKID_pt = lsim(sysERAOKID,u,tempo);

    % Figura - Resposta ao Impulso dos ROMs a partir do DS comparado com CFD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run('figures/figures_DSresp')
    close all
end