clearvars
close all
clc
addpath('figures')
addpath('functions')

%% Resposta ao Impulso a partir do Unit Sample (US)

% Dados do Fort.14 (US)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID1 = fopen([fileparts(pwd),filesep,'0_data/USplunge_Fort.14']);
fileID2 = fopen([fileparts(pwd),filesep,'0_data/USpitch_Fort.14']);
C1 = textscan(fileID1,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
C2 = textscan(fileID2,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID1);
fclose(fileID2);
[~,ht_pg,~,~,~,~,cl_pg,cm_pg,~] = deal(C1{:,1:9});
[t,~,aoa_pt,~,~,~,cl_pt,cm_pt,~] = deal(C2{:,1:9});
N = size(t,1);
for i=2:N
    yUS(:,:,i-1) = [-cl_pg(i) cl_pt(i); -cm_pg(i) cm_pt(i)];
end
yUS(:,:,N) = yUS(:,:,end);

% Salvar Resposta ao Impulso (US)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename =  [fileparts(pwd),filesep,'2_pipeline/yUS.mat'];
save(filename,'yUS');

%% Resposta ao Impulso a partir do Discrete Step (DS)

% Amplitudes Máximas em Plunge e Pitch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amplitude = ["low"; "med"; "high"];
for i=1:length(amplitude)
    amp = amplitude(i,:);
    switch amp
        case 'low'
            amp_plunge = '1e-6';
            amp_pitch = '1e-4';
            amp_pg_lowAmp = str2double(amp_plunge);
            amp_pt_lowAmp = str2double(amp_pitch);
        case 'med'
            amp_plunge = '500e-6';
            amp_pitch = '500e-4';
            amp_pg_medAmp = str2double(amp_plunge);
            amp_pt_medAmp = str2double(amp_pitch);
        case 'high'
            amp_plunge = '1000e-6';
            amp_pitch = '1000e-4';
            amp_pg_highAmp = str2double(amp_plunge);
            amp_pt_highAmp = str2double(amp_pitch);
    end
    amp_pg = str2double(amp_plunge);
    amp_pt = str2double(amp_pitch);

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
    N = length(t);

    % Convenção de Sinal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcl_pg = -dcl_pg;
    dcm_pg = -dcm_pg;

    % Resposta ao Impulso (DS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yDS = zeros(2,2,N);
    for i=2:N
        yDS(:,:,i-1) = [dcl_pg(i) dcl_pt(i); dcm_pg(i) dcm_pt(i)];
    end
    yDS(:,:,N) = yDS(:,:,end);

    % Salvar Resposta ao Impulso (DS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch amp
        case 'low'
            yDS_lowAmp = yDS;
            filename = '2_pipeline/ROMid/yDS_lowAmp.mat';
            save([fileparts(pwd),filesep,filename],'yDS_lowAmp');
        case 'med'
            yDS_medAmp = yDS;
            filename = '2_pipeline/ROMid/yDS_medAmp.mat';
            save([fileparts(pwd),filesep,filename],'yDS_medAmp');
        case 'high'
            yDS_highAmp = yDS;
            filename = '2_pipeline/ROMid/yDS_highAmp.mat';
            save([fileparts(pwd),filesep,filename],'yDS_highAmp');
    end

    % Figura - Resposta Temporal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    run('./figures/figures_yDS_time');
    close all
end