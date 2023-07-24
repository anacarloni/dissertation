% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars -except caso casos Np Npmax
format long
close all
addpath('functions');

% Caso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caso = input('Caso [01-28]: ','s');

% Quantidade de Polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Np = input('Quantidade de polos [1-6]: ');

% Aquisição de Dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = sprintf('2_pipeline/output_spectral/Case%s_fit_in_a.dat',caso);
filename2 = sprintf('2_pipeline/output_spectral/Case%s_fit_in_h.dat',caso);
fileID1 = fopen([fileparts(pwd),filesep,filename1]);
fileID2 = fopen([fileparts(pwd),filesep,filename2]);
C1 = textscan(fileID1,'%f %f %f %f %f');
C2 = textscan(fileID2,'%f %f %f %f %f');
fclose(fileID1);
fclose(fileID2);
[~,recla,imcla,recma,imcma] = deal(C1{:,1:5});
[k,reclh,imclh,recmh,imcmh] = deal(C2{:,1:5});
N = size(k,1);
kk = linspace(0,max(k),1001)';

% Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = struct('Clh',zeros(N,1),'Cla',zeros(N,1),'Cmh',zeros(N,1),'Cma',zeros(N,1));
M = struct('Clh',zeros(N,1),'Cla',zeros(N,1),'Cmh',zeros(N,1),'Cma',zeros(N,1));

% Função de Transferência
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H.Clh = complex(reclh,imclh);
H.Cla = complex(recla,imcla);
H.Cmh = complex(recmh,imcmh);
H.Cma = complex(recma,imcma);

% Módulo Quadrático
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M.Clh = abs(H.Clh).^2;
M.Cla = abs(H.Cla).^2;
M.Cmh = abs(H.Cmh).^2;
M.Cma = abs(H.Cma).^2;

for i = 1:N
    if(M.Clh(i) < 1)
        M.Clh(i) = 1;
    end
    if(M.Cla(i) < 1)
        M.Cla(i) = 1;
    end
    if(M.Cmh(i) < 1)
        M.Cmh(i) = 1;
    end
    if(M.Cma(i) < 1)
        M.Cma(i) = 1;
    end
end

% Conjunto de Polos Iniciais Igualmente Espaçados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lbd = 0.25;
ubd = 1.2;
Bini = linspace(lbd,ubd,Np);

% Coeficientes Lineares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,Fc,kFc,C] = fiterror_1stForm(Bini,H,M,k,kk);
Bf = Bini;

% Pontos do Polinômio para Comparação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLhfit = Fc*C.Clh;
CLafit = Fc*C.Cla;
CMhfit = Fc*C.Cmh;
CMafit = Fc*C.Cma;
kCLhfit = kFc*C.Clh;
kCLafit = kFc*C.Cla;
kCMhfit = kFc*C.Cmh;
kCMafit = kFc*C.Cma;

% Partes Reais e Imaginárias dos Coeficientes Lineares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reCclh = real(C.Clh).';
imCclh = imag(C.Clh).';
reCcmh = real(C.Cmh).';
imCcmh = imag(C.Cmh).';
reCcla = real(C.Cla).';
imCcla = imag(C.Cla).';
reCcma = real(C.Cma).';
imCcma = imag(C.Cma).';

% Comparação dos Resultados - Função da Frequência Reduzida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('figures');
dir = [fileparts(pwd),filesep,'3_output/unoptimized/' num2str(Np) '_poles\'];
figures_tf

% Comparação dos Resultados - Plano Complexo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kflutter = 0.1826;
kref = kflutter/2;
id = find(k>=kref,1);
kref = kflutter;
id2 = find(k>=kref,1);
switch caso
    case {'04'}
        aux=20;
    case {'05'}
        aux=10;
    case {'13','18'}
        aux = 8;
    otherwise
        aux = 16;
end
% figures_complex_tf

% Comparação de Resultados Pontuais - Plano Complexo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures_complex_fit_tf

% Norma L2 entre RFA e SIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L2norm

% Arquivos de Saída
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1 = Bf;
out2 = [reCclh', reCcmh', reCcla', reCcma'];
out3 = [k real(CLafit) imag(CLafit) real(CMafit) imag(CMafit)];
out4 = [k real(CLhfit) imag(CLhfit) real(CMhfit) imag(CMhfit)];
out5 = [clhL2_WF claL2_WF cmhL2_WF cmaL2_WF; clhL2_RFA claL2_RFA cmhL2_RFA cmaL2_RFA];

fout1 = fopen(sprintf('Case%s_poles.dat',caso),'w');
fout2 = fopen(sprintf('Case%s_coef.dat',caso),'w');
fout3 = fopen(sprintf('Case%s_fit_out_a.dat',caso),'w');
fout4 = fopen(sprintf('Case%s_fit_out_h.dat',caso),'w');
fout5 = fopen(sprintf('Case%s_L2norm.dat',caso),'w');

fprintf(fout1,'%11.7E\n',out1);

fprintf(fout2,'%11.7E %11.7E %11.7E %11.7E\n',out2');

fprintf(fout3,'VARIABLES = "k" "Re(Cl)" "Im(Cl)" "Re(Cm)" "Im(Cm)"\n');
fprintf(fout3,'ZONE\n');
fprintf(fout3,'I=%4d, J=1, K=1, F=BLOCK\n',N);
fprintf(fout3,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out3');

fprintf(fout4,'VARIABLES = "k" "Re(Cl)" "Im(Cl)" "Re(Cm)" "Im(Cm)"\n');
fprintf(fout4,'ZONE\n');
fprintf(fout4,'I=%4d, J=1, K=1, F=BLOCK\n',N);
fprintf(fout4,'%11.7E %11.7E %11.7E %11.7E %11.7E\n',out4');

fprintf(fout5,'VARIABLES = "Clh" "Cla" "Cmh" "Cma"\n');
fprintf(fout5,'%11.7E %11.7E %11.7E %11.7E\n',out5');

fclose(fout1);
fclose(fout2);
fclose(fout3);
fclose(fout4);
fclose(fout5);

filename = [fileparts(pwd),filesep,...
    sprintf('2_pipeline/unoptimized/%s_poles',num2str(Np))];
movefile('*.dat',filename)