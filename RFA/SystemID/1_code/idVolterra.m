%% Identificar ROM a partir da Resposta ao Impulso (DS) usando Volterra

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

    % Dados do Fort.14 (DS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    % Convenção de Sinal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dcl_pg = -dcl_pg;
    dcm_pg = -dcm_pg;

    % Variação dos Sinais de Entrada
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dht = ht;
    daoa = aoa;
    for j = N:-1:2
        dht(j,1) = ht(j,1)-ht(j-1,1);
        daoa(j,1) = aoa(j,1)-aoa(j-1,1);
    end
    dcl_pg = [0; dcl_pg];
    dcm_pg = [0; dcm_pg];
    dcl_pt = [0; dcl_pt];
    dcm_pt = [0; dcm_pt];

    % Normalização
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dht = dht./amp_pg;
    daoa = daoa./amp_pt;
    dcl_pg = dcl_pg./amp_pg;
    dcm_pg = dcm_pg./amp_pg;
    dcl_pt = dcl_pt./amp_pt;
    dcm_pt = dcm_pt./amp_pt;

    % Sinais de Entrada e Saída
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u1 = ht(3:1e3);
    du1 = dht(3:1e3);
    u2 = aoa(3:1e3);
    du2 = daoa(3:1e3);
    s11 = dcl_pg(3:1e3);
    s12 = dcm_pg(3:1e3);
    s21 = dcl_pt(3:1e3);
    s22 = dcm_pt(3:1e3);

    % Resposta a Sinal de Entrada Arbitrário com Primeiro Kernel de Volterra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y11 = (u1(1)*s11)' + conv(s11,du1);
    y12 = (u1(1)*s12)' + conv(s12,du1);
    y21 = (u2(1)*s21)' + conv(s21,du2);
    y22 = (u2(1)*s22)' + conv(s22,du2);

    % Resposta a Sinal de Entrada Arbitrário para Cada Amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch amp
        case 'low'
            y11_low = y11;
            y12_low = y12;
            y21_low = y21;
            y22_low = y22;
            filename = '2_pipeline/idVolterra/Volterra_ylow.mat';
            save([fileparts(pwd),filesep,filename],'y11_low','y12_low','y21_low','y22_low');
            filename = '2_pipeline/idVolterra/Volterra_sLow.mat';
            save([fileparts(pwd),filesep,filename],'s11','s12','s21','s22');
        case 'med'
            y11_med = y11;
            y12_med = y12;
            y21_med = y21;
            y22_med = y22;
            filename = '2_pipeline/idVolterra/Volterra_ymed.mat';
            save([fileparts(pwd),filesep,filename],'y11_med','y12_med','y21_med','y22_med');
        case 'high'
            y11_high = y11;
            y12_high = y12;
            y21_high = y21;
            y22_high = y22;
            filename = '2_pipeline/idVolterra/Volterra_yhigh.mat';
            save([fileparts(pwd),filesep,filename],'y11_high','y12_high','y21_high','y22_high');
    end

    % Resposta a Sinal de Entrada Arbitrário para Baixa Amplitude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch amp
        case 'low'
            yVolt = zeros(2,2,N);
            for i=1:N
                yVolt(:,:,i) = [y11_low(i)*amp_pg y21_low(i)*amp_pt;
                    y12_low(i)*amp_pg y22_low(i)*amp_pt];
            end
            filename = '2_pipeline/idVolterra/Volterra_yVolt.mat';
            save([fileparts(pwd),filesep,filename],'yVolt');
    end
end

%% Resposta ao Impulso com Primeiro Kernel de Volterra

% Abrir Resposta a Sinal de Entrada Arbitrário para Cada Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_sLow.mat'])
load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_ylow.mat'])
load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_ymed.mat'])
load([fileparts(pwd),filesep,'2_pipeline/idVolterra/Volterra_yhigh.mat'])

% Figuras da Resposta ao Impulso (DS) e Primeiro Kernel de Volterra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('figures/figures_impulseDS_Volt')