% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars -except inputsignal nblocks caso casos
format long e
close all
addpath('functions');
addpath('figures');

% Número de Mach do Escoamento Não Perturbado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 0.8;

% Frequência Reduzida Selecionada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = 0.2;

% Aquisição de Dados do Arquivo Fort.14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputsignal = input('Sinal de entrada [WF1-WF7]: ','s');
caso = input('Caso [01-28]: ','s');
filename = sprintf('0_data/Input%s_Fort.14',inputsignal);
fileID = fopen([fileparts(pwd),filesep,filename]);
C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID);
[t,ht,aoa,dcl,cd,dcm,cl,cm,cmex] = deal(C{:,1:9});
N = size(t,1);
Ts = t(2)-t(1);
Fs = 1/Ts;

% Detectando Pulso(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hflag = max(abs(ht));
aflag = max(abs(aoa));
if (hflag == 0)
    hflag = 0;
else
    hflag = 1;
end
if (aflag == 0)
    aflag = 0;
else
    aflag = 1;
end

% Convenção e Unidades
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ht = -ht;
aoa = pi/180*aoa;

% Derivada dos Sinais
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dht = ht;
daoa = aoa;
for j = N:-1:2
    dht(j,1) = ht(j,1)-ht(j-1,1);
    daoa(j,1) = aoa(j,1)-aoa(j-1,1);
end
dcl(1) = 0;
dcm(1) = 0;

% Número de Blocos no Sinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nblocks = input('Número de blocos: ');
nfft = floor(N/nblocks);

% Extraindo o Valor Inicial nos Blocos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
htz = ht;
aoaz = aoa;
clz = cl;
cmz = cm;
for jj = 1:nblocks
    jk1 = (jj-1)*nfft+1;
    jk2 = jj*nfft;
    for jk = jk1:jk2
        htz(jk,1) = ht(jk,1)-ht(jk1,1);
        aoaz(jk,1) = aoa(jk,1)-aoa(jk1,1);
        clz(jk,1) = cl(jk,1)-cl(jk1,1);
        cmz(jk,1) = cm(jk,1)-cm(jk1,1);
    end
end
htz(N,1) = htz(N-1,1);
aoaz(N,1) = aoaz(N-1,1);
clz(N,1) = clz(N-1,1);
cmz(N,1) = cmz(N-1,1);

% Entradas e Saídas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = dht;
u2 = daoa;
y1 = dcl;
y2 = dcm;

% Correlação entre Sinais e Derivadas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produto interno
prod_p = dot(ht,aoa);
prod_dp = dot(dht,daoa);

% Integral
int_p = trapz(ht.*aoa);
int_dp = trapz(dht.*daoa);

% Cross power spectral density
% figure
% cpsd(ht,aoa)
% title('CPSD (h, AoA)')
% figure
% cpsd(u1,u2)
% title('CPSD (dh, dAoA)')

% Densidade Espectral de Potência (PSD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choose_case
novlp = novlp*nfft;
nfft = nfft*nsample;
% window = input('Janela [rectwin,hanning]: ','s');
switch window
    case 'rectwin'
        wind = rectwin(nfft);
    case 'hanning'
        wind = hanning(nfft);
end
[G11,~] = tfestimate(u1,y1,wind,novlp,nfft,Fs);
[G21,~] = tfestimate(u1,y2,wind,novlp,nfft,Fs);
[G12,~] = tfestimate(u2,y1,wind,novlp,nfft,Fs);
[G22,f] = tfestimate(u2,y2,wind,novlp,nfft,Fs);

% Frequência Angular e Reduzida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 2*pi*f;
k = w./(2*M);

% Aproximação para a Frequência Selecionada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dk = abs(k-ks);
dkmax = max(dk);
[~,iks] = min(dk);
ktemp1 = k(iks);
G11temp1 = G11(iks);
G12temp1 = G12(iks);
G21temp1 = G21(iks);
G22temp1 = G22(iks);

[~,iks] = max(iks);
ktemp2 = k(iks);
G11temp2 = G11(iks);
G12temp2 = G12(iks);
G21temp2 = G21(iks);
G22temp2 = G22(iks);

G11s = (ks-ktemp1)*(G11temp2-G11temp1)/(ktemp2-ktemp1)+G11temp1;
G12s = (ks-ktemp1)*(G12temp2-G12temp1)/(ktemp2-ktemp1)+G12temp1;
G21s = (ks-ktemp1)*(G21temp2-G21temp1)/(ktemp2-ktemp1)+G21temp1;
G22s = (ks-ktemp1)*(G22temp2-G22temp1)/(ktemp2-ktemp1)+G22temp1;

% Módulos, Fases, Partes Reais e Imaginárias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mag11 = abs(G11);
phi11 = angle(G11);
Gr11 = real(G11);
Gi11 = imag(G11);

mag21 = abs(G21);
phi21 = angle(G21);
Gr21 = real(G21);
Gi21 = imag(G21);

mag12 = abs(G12);
phi12 = angle(G12);
Gr12 = real(G12);
Gi12 = imag(G12);

mag22 = abs(G22);
phi22 = angle(G22);
Gr22 = real(G22);
Gi22 = imag(G22);

mag11s = abs(G11s);
phi11s = angle(G11s);
Gr11s = real(G11s);
Gi11s = imag(G11s);

mag21s = abs(G21s);
phi21s = angle(G21s);
Gr21s = real(G21s);
Gi21s = imag(G21s);

mag12s = abs(G12s);
phi12s = angle(G12s);
Gr12s = real(G12s);
Gi12s = imag(G12s);

mag22s = abs(G22s);
phi22s = angle(G22s);
Gr22s = real(G22s);
Gi22s = imag(G22s);

% Diagramas de Bode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(2,1,1);
% plot(k,log(mag11))
% ylabel('Magnitude (dB)')
% title('G_C_l_,_h')
% grid on
% subplot(2,1,2);
% plot(k,phi11)
% ylabel('Phase (rad)')
% xlabel('\kappa')
% grid on
%
% figure
% subplot(2,1,1);
% plot(k,log(mag21))
% ylabel('Magnitude (dB)')
% title('G_C_m_,_h')
% grid on
% subplot(2,1,2);
% plot(k,phi21)
% ylabel('Phase (rad)')
% xlabel('\kappa')
% grid on
%
% figure
% subplot(2,1,1);
% plot(k,log(mag12))
% ylabel('Magnitude (dB)')
% title('G_C_l_,_\alpha')
% grid on
% subplot(2,1,2);
% plot(k,phi12)
% ylabel('Phase (rad)')
% xlabel('\kappa')
% grid on
%
% figure
% subplot(2,1,1);
% plot(k,log(mag22))
% ylabel('Magnitude (dB)')
% title('G_C_m_,_\alpha')
% grid on
% subplot(2,1,2);
% plot(k,phi22)
% ylabel('Phase (rad)')
% xlabel('\kappa')
% grid on

% Partes Reais e Imaginárias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(k,Gr11,'b',k,Gi11,'g')
% % xlim([0 1])
% % ylim([-.2 7])
% ylabel('G_C_l_,_h')
% xlabel('\kappa')
% legend('Re','Im')
% grid on
%
% figure
% plot(k,Gr12,'b',k,Gi12,'g')
% % xlim([0 1])
% % ylim([-4.5 14.5])
% ylabel('G_C_l_,_\alpha')
% xlabel('\kappa')
% legend('Re','Im')
% grid on
%
% figure
% plot(k,Gr21,'b',k,Gi21,'g')
% % xlim([0 1])
% % ylim([-1.3 .9])
% ylabel('G_C_m_,_h')
% xlabel('\kappa')
% legend('Re','Im')
% grid on
%
% figure
% plot(k,Gr22,'b',k,Gi22,'g')
% % xlim([0 1])
% % ylim([-1.5 .6])
% ylabel('G_C_m_,_\alpha')
% xlabel('\kappa')
% legend('Re','Im')
% grid on

% Dados do Caso 00 (SIMO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(caso,'04') || strcmp(caso,'11') || strcmp(caso,'12') || strcmp(caso,'16') || strcmp(caso,'17'))
        overlap = '1';
elseif (strcmp(caso,'05') || strcmp(caso,'13') || strcmp(caso,'18'))
        overlap = '0';
end
aux = sprintf('0_data/fit_in_h_c00_case%s.dat',caso);
fileID = fopen([fileparts(pwd),filesep,aux]);
C = textscan(fileID,'%f %f %f %f %f');
fclose(fileID);
[~,Gr01_c00,Gi01_c00,Gr03_c00,Gi03_c00] = deal(C{:,1:5});

aux = sprintf('0_data/fit_in_a_c00_case%s.dat',caso);
fileID = fopen([fileparts(pwd),filesep,aux]);
C = textscan(fileID,'%f %f %f %f %f');
fclose(fileID);
[kk_c00,Gr02_c00,Gi02_c00,Gr04_c00,Gi04_c00] = deal(C{:,1:5});

% Sobreposição das Partes Reais e Imaginárias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures_intpol

% Vetores de Saída
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(k,1);
klim = 3;
for j = 1:N
    if k(j) <= klim
        kk(j,1) = k(j,1);
        Gr01(j,1) = Gr11(j,1);
        Gr02(j,1) = Gr12(j,1);
        Gr03(j,1) = Gr21(j,1);
        Gr04(j,1) = Gr22(j,1);
        Gi01(j,1) = Gi11(j,1);
        Gi02(j,1) = Gi12(j,1);
        Gi03(j,1) = Gi21(j,1);
        Gi04(j,1) = Gi22(j,1);
        mag01(j,1) = mag11(j,1);
        mag02(j,1) = mag12(j,1);
        mag03(j,1) = mag21(j,1);
        mag04(j,1) = mag22(j,1);
        phi01(j,1) = phi11(j,1);
        phi02(j,1) = phi12(j,1);
        phi03(j,1) = phi21(j,1);
        phi04(j,1) = phi22(j,1);
    end
end
N = size(kk,1);

% Arquivos de Saída
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1 = [kk Gr01 Gi01 Gr02 Gi02 Gr03 Gi03 Gr04 Gi04];
out2 = [kk mag01 phi01 mag02 phi02 mag03 phi03 mag04 phi04];
if (hflag == 1)
    out3 = [kk Gr01 Gi01 Gr03 Gi03];
end
if (aflag == 1)
    out4 = [kk Gr02 Gi02 Gr04 Gi04];
end

out1 = out1';
out2 = out2';
if (hflag == 1)
    out3 = out3';
end
if (aflag == 1)
    out4 = out4';
end

if (hflag == 0)
    fout1 = fopen(sprintf('Case%s_FFT1_a.dat',caso),'w');
    fout2 = fopen(sprintf('Case%s_FFT2_a.dat',caso),'w');
elseif (aflag == 0)
    fout1 = fopen(sprintf('Case%s_FFT1_h.dat',caso),'w');
    fout2 = fopen(sprintf('Case%s_FFT2_h.dat',caso),'w');
else
    fout1 = fopen(sprintf('Case%s_FFT1.dat',caso),'w');
    fout2 = fopen(sprintf('Case%s_FFT2.dat',caso),'w');
end
if (hflag == 1)
    fout3 = fopen(sprintf('Case%s_fit_in_h.dat',caso),'w');
end
if (aflag == 1)
    fout4 = fopen(sprintf('Case%s_fit_in_a.dat',caso),'w');
end

fprintf(fout1,'VARIABLES = "k" "Re(CL)_h" "Im(CL)_h" "Re(CL)_a" "Im(CL)_a" "Re(CM)_h" "Im(CM)_h" "Re(CM)_a" "Im(CM)_a"\n');
fprintf(fout1,'ZONE\n');
fprintf(fout1,'I=%4d, J=1, K=1, F=BLOCK\n',N);
fprintf(fout1,'%11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E\n',out1);

fprintf(fout2,'VARIABLES = "k" "Abs(CL)_h" "Phase(CL)_h" "Abs(CL)_a" "Phase(CL)_a" "Abs(CM)_h" "Phase(CM)_h" "Abs(CM)_a" "Phase(CM)_a"\n');
fprintf(fout2,'ZONE\n');
fprintf(fout2,'I=%4d, J=1, K=1, F=BLOCK\n',N);
fprintf(fout2,'%11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E %11.7E\n',out2);

if (hflag == 1)
    fprintf(fout3,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out3);
end
if (aflag == 1)
    fprintf(fout4,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out4);
end

fclose(fout1);
fclose(fout2);
if (hflag == 1)
    fclose(fout3);
end
if (aflag == 1)
    fclose(fout4);
end

filename = [fileparts(pwd),filesep,'2_pipeline/output_spectral'];
movefile('*.dat',filename)

% Figuras dos Sinais de Entrada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = sprintf('0_data/meshmove_%s.dat',inputsignal);
fileID = fopen([fileparts(pwd),filesep,filename]);
C = textscan(fileID,'%f %f');
fclose(fileID);
[h_inp,aoa_inp] = deal(C{:,1:2});
N = size(t,1);

% figure
% subplot(2,1,1)
% plot(linspace(1,N,1000),h_inp(1:100:end)'./(abs(max(h_inp))),'b')
% ylabel('h/|h_{max}|')
% xlabel('Time Step')
% ylim([-1.5 1.5])
% xlim([0 N])
% grid on
% set(gca,'FontSize',17);
% 
% subplot(2,1,2)
% plot(linspace(1,N,1000),aoa_inp(1:100:end)'./(abs(max(aoa_inp))),'b')
% ylabel('\alpha/|\alpha_{max}|')
% xlabel('Time Step')
% ylim([-1.5 1.5])
% xlim([0 N])
% grid on
% set(gca,'FontSize',17);
%
% figname = sprintf('input_%s',inputsignal);
% hgsave([figname,'.fig'])
% print(figname,'-dpdf')
% filepath = [fileparts(pwd),filesep,'3_output/output_spectral/'];
% movefile(strcat(figname,'.pdf'),strcat(filepath,figname,'.pdf'))
% movefile(strcat(figname,'.fig'),strcat(filepath,figname,'.fig'))

% Figuras block/sample/overlap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures_systemvariables