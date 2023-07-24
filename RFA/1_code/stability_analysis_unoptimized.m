% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars -except caso casos Np Npmax
format long
close all
addpath('figures')

% Parâmetros do Sistema Aeroelástico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xa = 1.8;   % distance from the EA to the CM, normalized by semichord
xea = -1/4; % EA position used to calculate aerodynamic coefficients
ah = -2;    % distance from midchord to the EA, normalized by semichord
ra = 1.865; % airfoil dimensionless radius of gyration about the EA
wa = 100;   % uncoupled natural circular frequency of the pitch DOF, rad/s
wh = 100;   % uncoupled natural circular frequency of the plunge DOF, rad/s
wr = wa;    % reference circular frequency, rad/s
mu = 60;    % mass ratio, mu = m/(pi*rho_inf*b^2)

% Pressão Dinâmica Característica, Q*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qmin = 0;
Qmax = 1;
NQ = 100;
Q = linspace(Qmin,Qmax,NQ);

% Caso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caso = input('Caso [01-28]: ','s');

% Quantidade de Polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Np = input('Quantidade de polos [1-6]: ');

% Polos Não Otimizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form = 'unoptimized';

% Polos e Coeficientes Lineares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = sprintf('2_pipeline/%s/%d_poles/Case%s_poles.dat',...
    form,Np,caso);
filename2 = sprintf('2_pipeline/%s/%d_poles/Case%s_coef.dat',...
    form,Np,caso);
fileID1 = fopen([fileparts(pwd),filesep,filename1]);
fileID2 = fopen([fileparts(pwd),filesep,filename2]);
C1 = textscan(fileID1,'%f');
C2 = textscan(fileID2,'%f %f %f %f');
fclose(fileID1);
fclose(fileID2);
B = deal(C1{:,1});
[Cclh,Ccmh,Ccla,Ccma] = deal(C2{:,1:4});

% Corrigindo os Coeficientes (Posição do Eixo Elástico)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ccmh = -Ccmh;
arm = ah/2 - xea;
Ccma = Ccma + arm*(Ccla - Ccmh - arm*Cclh);
Ccla = Ccla - arm*Cclh;
Ccmh = Ccmh + arm*Cclh;

% Matrizes Aerodinâmica, de Massa e de Rigidez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(2,2,3+Np);
for i = 1:3+Np
    A(:,:,i) = [-Cclh(i)/2 -Ccla(i) ; Ccmh(i) 2*Ccma(i)];
end
M = [1 xa ; xa (ra^2)];
M2 = M - 1/(pi*mu)*A(:,:,3);
Minv = M\eye(size(M));
M2inv = M2\eye(size(M2));
K = [((wh/wr)^2) 0 ; 0 ((ra*wa/wr)^2)];

% Alocação de Memória
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flutter_i_U = 0;
U = zeros(1,NQ);
kappa = zeros(1,NQ);
sizeD = 2*(2+Np);
flutter_sol(1:NQ) = struct('eigval',zeros(sizeD,1),...
    'eigvec',zeros(sizeD,sizeD),...
    'D',zeros(sizeD,sizeD),...
    'wn',zeros(sizeD,1),...
    'fn',zeros(sizeD,1),...
    'zeta',zeros(sizeD,1));

% Matriz Dinâmica ou de Estabilidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:NQ
    U(i) = sqrt(Q(i)*mu);
    K2 = K - (U(i)^2)/(pi*mu)*A(:,:,1);
    C = -(U(i))/(pi*mu)*A(:,:,2);

    for j = 1:Np
        A2(:,:,j) = (U(i)^3)/(pi*mu)*A(:,:,(3+j));
    end
    if (Q(i) == 0)
        D(1:2,1:4) = [-Minv*C -Minv*K2];
    else
        D(1:2,1:4) = [-M2inv*C -M2inv*K2];
    end
    D(3:4,1:2) = eye(2);
    for j = 1:Np
        l = 5 + (j-1)*2;
        D(1:2,l:(l+1)) = M2inv*A2(:,:,j);
        D(l:(l+1),3:4) = eye(2);
        D(l:(l+1),l:(l+1)) = -U(i)*B(j)*eye(2);
    end
    kappa(i) = cond(D);

    % Problema de Autovalor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Vt,Vl] = eig(D);
    R(:,:,i) = Vt;
    S(:,i) = diag(Vl);
    eigvec = R(:,:,i);
    eigval = S(:,i);
    [eigval,i_sort] = sort(eigval);
    eigvec = eigvec(:,i_sort);
    [wn,zeta] = damp(D);

    flutter_sol(i).D = D;
    flutter_sol(i).eigvec = eigvec;
    flutter_sol(i).eigval = eigval;
    flutter_sol(i).wn = wn;
    flutter_sol(i).fn = wn/(2*pi);
    flutter_sol(i).zeta = zeta;

    if any(real(S(:,i))>1e-8) && flutter_i_U==0
        flutter_i_U = i;
    end
end

% Velocidade Crítica de Flutter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flutter_i_U>0
    min_flutter_speed = U(flutter_i_U);
end
disp('Velocidade característica de flutter:');
disp(min_flutter_speed);

% Autovalor e Autovetor de Flutter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('min_flutter_speed','var')
    flutter_id = find(flutter_sol(U==min_flutter_speed).eigval>0,1);
    eigval_flutter = flutter_sol(U==min_flutter_speed).eigval(flutter_id);
    eigvec_flutter = flutter_sol(U==min_flutter_speed).eigvec(:,flutter_id);
    disp('Autovalor de flutter:');
    disp(eigval_flutter);

    % Velocidade Crítica de Flutter - Interpolação quadrática
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flutter_i_U==1
    eigval_ant = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
    eigval_pos = flutter_sol(U==U(flutter_i_U+2)).eigval(flutter_id);
    x = [real(eigval_flutter), real(eigval_ant), real(eigval_pos)];
    y = [U(flutter_i_U), U(flutter_i_U+1), U(flutter_i_U+2)];
    coefficients = polyfit(x, y, 2);
    Uf_intquad = polyval(coefficients, 0);
    else
    eigval_ant = flutter_sol(U==U(flutter_i_U-1)).eigval(flutter_id);
    eigval_pos = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
    x = [real(eigval_ant), real(eigval_flutter), real(eigval_pos)];
    y = [U(flutter_i_U-1), U(flutter_i_U), U(flutter_i_U+1)];
    coefficients = polyfit(x, y, 2);
    Uf_intquad = polyval(coefficients, 0);
    end
    disp('Velocidade característica de flutter (interpolação quadrática):');
    disp(Uf_intquad);

    % Coeficientes de Amortecimento e Frequências Naturais
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    id = find(imag(flutter_sol(flutter_i_U).eigval)>4.5);
    id = id(1);
    if flutter_i_U==1
        eigval_ant = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
        eigval_pos = flutter_sol(U==U(flutter_i_U+2)).eigval(flutter_id);
        x = [imag(eigval_flutter), imag(eigval_ant), imag(eigval_pos)];
        y = [U(flutter_i_U), U(flutter_i_U+1), U(flutter_i_U+2)];
        coefficients = polyfit(y, x, 2);
        wn1 = polyval(coefficients, Uf_intquad);
        x = [real(eigval_flutter), real(eigval_ant), real(eigval_pos)];
        coefficients = polyfit(y, x, 2);
        sigma1 = polyval(coefficients, Uf_intquad);

        eigval_ant = flutter_sol(flutter_i_U+1).eigval(id);
        eigval_ft = flutter_sol(flutter_i_U).eigval(id);
        eigval_pos = flutter_sol(flutter_i_U+2).eigval(id);
        x = [imag(eigval_ft), imag(eigval_ant), imag(eigval_pos)];
        y = [U(flutter_i_U), U(flutter_i_U+1), U(flutter_i_U+2)];
        coefficients = polyfit(y, x, 2);
        wn2 = polyval(coefficients, Uf_intquad);
        x = [real(eigval_ft), real(eigval_ant), real(eigval_pos)];
        coefficients = polyfit(y, x, 2);
        sigma2 = polyval(coefficients, Uf_intquad);
    else
        eigval_ant = flutter_sol(U==U(flutter_i_U-1)).eigval(flutter_id);
        eigval_pos = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
        x = [imag(eigval_ant), imag(eigval_flutter), imag(eigval_pos)];
        y = [U(flutter_i_U-1), U(flutter_i_U), U(flutter_i_U+1)];
        coefficients = polyfit(y, x, 2);
        wn1 = polyval(coefficients, Uf_intquad);
        x = [real(eigval_ant), real(eigval_flutter), real(eigval_pos)];
        coefficients = polyfit(y, x, 2);
        sigma1 = polyval(coefficients, Uf_intquad);

        temp = find(imag(flutter_sol(flutter_i_U-1).eigval)>4.5);
        eigval_ant = flutter_sol(flutter_i_U-1).eigval(temp);
        eigval_ft = flutter_sol(flutter_i_U).eigval(id);
        eigval_pos = flutter_sol(flutter_i_U+1).eigval(id);
        x = [abs(imag(eigval_ant)), abs(imag(eigval_ft)), abs(imag(eigval_pos))];
        y = [U(flutter_i_U-1), U(flutter_i_U), U(flutter_i_U+1)];
        coefficients = polyfit(y, x, 2);
        wn2 = polyval(coefficients, Uf_intquad);
        x = [real(eigval_ant), real(eigval_ft), real(eigval_pos)];
        coefficients = polyfit(y, x, 2);
        sigma2 = polyval(coefficients, Uf_intquad);
    end
    sigma_wn = [sigma1 abs(wn1); sigma2 abs(wn2)];
    disp('Amortecimento e Frequência natural:');
    disp(sigma_wn);
end

% Separando os Autovetores/Valores Conjugados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sizeS,~] = size(S);
sizeS = sizeS/2;
for i = 1:sizeS
    R1(:,i,:) = R(:,2*i-1,:);
    R2(:,i,:) = R(:,2*i-0,:);
    S1(i,:) = S(2*i-1,:);
    S2(i,:) = S(2*i-0,:);
end

% Alinhando os Autovalores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:NQ
    R3 = R1(:,:,i);
    R4 = R2(:,:,i);
    S3 = S1(:,i);
    S4 = S2(:,i);
    dR = eye(sizeS);
    dS = eye(sizeS);
    dRS = eye(sizeS);
    for j = 1:sizeS
        for k = 1:sizeS
            dR(j,k) = sum(abs(R1(:,k,i)-R1(:,j,i-1)));
            dS(j,k) = abs(S1(k,i)-S1(j,i-1));
            dRS(j,k) = 1/(1/dR(j,k)^2+1/dS(j,k)^2);
        end
    end
    dRSmax = max(max(dRS))+1;
    for j = 1:sizeS
        [dRSmin,i1] = min(dRS);
        [dRSmin,i2] = min(dRSmin);
        i1 = i1(i2);
        dRS(i1,:) = dRSmax;
        dRS(:,i2) = dRSmax;
        i3(j,:) = [i1 i2];
    end
    i3 = sortrows(i3);
    for j = 1:sizeS
        R1(:,j,i) = R3(:,i3(j,2));
        R2(:,j,i) = R4(:,i3(j,2));
        S1(j,i) = S3(i3(j,2));
        S2(j,i) = S4(i3(j,2));
    end
end

% Isolando os Modos Estruturais
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear S3 S4
j = 0;
for i = 1:sizeS
    if (S1(i,1) ~= 0)
        j = j+1;
        S3(j,:) = S1(i,:);
        S4(j,:) = S2(i,:);
    end
end

% Resultados de Referência (Rausch, Batina and Yang, 1989)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_ref = [.2; .5; .8];
sigma_ref = [-.011 .0001 .017 -.068 -.148 -.223];
wn_ref = [.790 .913 1.022 5.353 5.349 5.317];

% Resultados do Caso 00 (Discrete Step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([fileparts(pwd),filesep,'0_data/S_c00.mat']);
aux = NQ/10;
rec00 = real(S_c00);
imc00 = imag(S_c00);

% Local das Raízes do Primeiro Modo Aeroelástico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = sprintf('3_output/stability_analysis/%s/%s_poles/',...
    form,num2str(Np));
dir = [fileparts(pwd),filesep,filename];
figures_rlocus_1stMode

% Local das Raízes do Segundo Modo Aeroelástico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figures_rlocus_2ndMode

% Arquivos de Saída
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep1 = real(S4(1,:));
imp1 = imag(S4(1,:));
rep2 = real(S4(2,:));
imp2 = imag(S4(2,:));

out = [Q; U; rep1; imp1; rep2; imp2];
fout = fopen(sprintf('Case%s_root.dat',caso),'w');
fprintf(fout,'VARIABLES = "Q" "U" "s1" "w1" "s2" "w2"\n');
fprintf(fout,'ZONE\n');
fprintf(fout,'I=%4d, J=1, K=1, F=POINT\n',NQ);
fprintf(fout,'%11.7E %11.7E %11.7E %11.7E %11.7E %11.7E\n',out);
fclose(fout);

if exist('min_flutter_speed','var')
    out = [Q(flutter_i_U); Uf_intquad; sigma_wn(1,1); sigma_wn(1,2);
        sigma_wn(2,1); sigma_wn(2,2)];
    fout = fopen(sprintf('Case%s_flutter.dat',caso),'w');
    fprintf(fout,'VARIABLES = "Q" "U" "s1" "w1" "s2" "w2" \n');
    fprintf(fout,'ZONE\n');
    fprintf(fout,'I=%4d, J=1, K=1, F=POINT\n',NQ);
    fprintf(fout,'%11.7E %11.7E %11.7E %11.7E %11.7E %11.7E\n',out');
    fclose(fout);
end

fout = fopen(sprintf('Case%s_full.dat',caso),'w');
for i = 1:sizeS
    rep1 = real(S1(i,:));
    imp1 = imag(S1(i,:));
    out = [rep1; imp1; U];
    fprintf(fout,'VARIABLES = "Real" "Imag" "U"\n');
    fprintf(fout,'ZONE\n');
    fprintf(fout,'I=%4d, J=1, K=1, F=POINT\n',NQ);
    fprintf(fout,'%11.7E %11.7E %11.7E\n',out);
end
fclose(fout);

filename = [fileparts(pwd),filesep,...
    sprintf('2_pipeline/stability_analysis/%s/%s_poles/',...
    form,num2str(Np))];
movefile('*.dat',filename)