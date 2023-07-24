% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
format long e
close all

% Número de Mach do Escoamento
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 0.8;

% Seleção de Frequência
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ks = 0.2;

% Aquisição de Dados do 'Fort.14'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputsignal = input('Sinal de entrada [WF1-WF7]: ','s');
signal = input('Plunge or pitch (h/a): ','s');
filename = sprintf('InputDS_%s_Fort.14',signal);
fileID = fopen(filename);
C = textscan(fileID,'%f %f %f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID);
[t,ht,aoa,dcl,cd,dcm,cl,cm,cmex] = deal(C{:,1:9});
N = size(t,1);
Ts = t(2)-t(1);
Fs = 1/Ts;

% Detectando o(s) Pulso(s)
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

% Extraíndo o valor Inicial nos Blocos
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

% Densidade Espectral de Potência
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = input('Caso [01-28]: ','s');
choose_case
novlp = novlp*nfft;
nfft = nfft*nsample;
window = input('Janela [rectwin,hanning]: ','s');
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

% Frequencia Angular e Reduzida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 2*pi*f;
k = w./(2*M);

% Aproximação para a Frequência Selecionada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dk = abs(k-ks);
[~,iks] = min(dk);
ktemp1 = k(iks);
G11temp1 = G11(iks);
G12temp1 = G12(iks);
G21temp1 = G21(iks);
G22temp1 = G22(iks);

[~,iks] = max(dk);
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

% Funções de Transferência
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure_intpol

% Arquivos de Saída
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (hflag == 1)
    out3 = [kk Gr01 Gi01 Gr03 Gi03];
end
if (aflag == 1)
    out4 = [kk Gr02 Gi02 Gr04 Gi04];
end

if (hflag == 1)
    out3 = out3';
end
if (aflag == 1)
    out4 = out4';
end

if (hflag == 1)
    fout3 = fopen(sprintf('fit_in_%s_c00_case%s.dat',signal,caso),'w');
end
if (aflag == 1)
    fout4 = fopen(sprintf('fit_in_%s_c00_case%s.dat',signal,caso),'w');
end

if (hflag == 1)
    fprintf(fout3,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out3);
end
if (aflag == 1)
    fprintf(fout4,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out4);
end

if (hflag == 1)
    fclose(fout3);
end
if (aflag == 1)
    fclose(fout4);
end
