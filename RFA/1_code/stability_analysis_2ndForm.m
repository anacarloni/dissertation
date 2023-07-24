% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars -except caso casos Np Npmax
format long
close all
addpath('functions')

% Caso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% caso = input('Caso [01-28]: ','s');

% Quantidade de Polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Np = input('Quantidade de polos [4-6]: ');

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
b = 1;      % semichord
ndof = 2;   % number of degrees of freedom, pitch and plunge

% Pressão Dinâmica Característica, Q*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qmin = 0;
Qmax = 1;
NQ = 100;
Q = linspace(Qmin,Qmax,NQ);

% Segunda Forma dos RFAs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
form = '2nd';

% Polos e Coeficientes Lineares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = sprintf('2_pipeline/optimized_%sForm/%d_poles/Case%s_poles.dat',...
    form,Np,caso);
filename2 = sprintf('2_pipeline/optimized_%sForm/%d_poles/Case%s_coef.dat',...
    form,Np,caso);
fileID1 = fopen([fileparts(pwd),filesep,filename1]);
fileID2 = fopen([fileparts(pwd),filesep,filename2]);
C1 = textscan(fileID1,'%f');
C2 = textscan(fileID2,'%f %f %f %f');
fclose(fileID1);
fclose(fileID2);
Bini = deal(C1{:,1});
[Cclh,Ccmh,Ccla,Ccma] = deal(C2{:,1:4});

% Polos com Multiplicidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,mult,p,index] = multiplicity2(Bini');
nb = Np;                        % número total de polos
N1 = p(1);                      %  N1 polos c/ multiplicidade 1
N2 = p(1)+p(2);                 % (N2-N1) polos c/ multiplicidade 2
N3 = p(1)+p(2)+p(3);            % (N3-N2) polos c/ multiplicidade 3
N4 = p(1)+p(2)+p(3)+p(4);       % (N4-N3) polos c/ multiplicidade 4
nb1 = length(find(mult==1));
nb2 = length(find(mult==2));
nb3 = length(find(mult==3));
nb4 = length(find(mult==4));

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
Uast = zeros(1,NQ);
kappa = zeros(1,NQ);
sizeD = 2*(2+nb);
D = zeros(sizeD);
A1tilde = zeros(2,2,nb1);
A2tilde = zeros(2,2,nb2);
A3tilde = zeros(2,2,nb3);
A4tilde = zeros(2,2,nb4);
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

    for j=1:nb1
        A1tilde(:,:,j) = (U(i)^3)/(pi*mu)*A(:,:,3+j); %q(i)*(U(i)/b)*
    end
    if nb1+1<=nb
        for j=nb1+1:nb1+nb2
            A2tilde(:,:,j-nb1) = (U(i)^4)/(pi*mu)*A(:,:,3+j);
        end
    end
    if nb1+nb2+1<=nb
        for j=nb1+nb2+1:nb1+nb2+nb3
            A3tilde(:,:,j-nb2-nb1) = (U(i)^5)/(pi*mu)*A(:,:,3+j);
        end
    end
    if nb1+nb2+nb3+1<=nb
        for j=nb1+nb2+nb3+1:nb1+nb2+nb3+nb4
            A4tilde(:,:,j-nb1-nb2-nb3) = (U(i)^6)/(pi*mu)*A(:,:,3+j);
        end
    end

    if (Q(i) == 0)
        D(1:ndof,1:2*ndof) = [-Minv*C -Minv*K2];
        D(1:ndof,2*ndof+1:2*(ndof+nb1)) = Minv*reshape(A1tilde,[ndof,2*nb1]);
        if nb1+1<=nb
            D(1:ndof,2*(ndof+nb1)+1:2*(ndof+nb1+nb2)) = Minv*reshape(A2tilde,[ndof,2*nb2]);
        end
        if nb1+nb2+1<=nb
            D(1:ndof,2*(ndof+nb1+nb2)+1:2*(ndof+nb1+nb2+nb3)) = Minv*reshape(A3tilde,[ndof,2*nb3]);
        end
        if nb1+nb2+nb3+1<=nb
            D(1:ndof,2*(ndof+nb1+nb2+nb3)+1:2*(ndof+nb1+nb2+nb3+nb4)) = Minv*reshape(A4tilde,[ndof,2*nb4]);
        end
    else
        D(1:ndof,1:2*ndof) = [-M2inv*C -M2inv*K2];
        D(1:ndof,2*ndof+1:2*(ndof+nb1)) = M2inv*reshape(A1tilde,[ndof,2*nb1]);
        if nb1+1<=nb
            D(1:ndof,2*(ndof+nb1)+1:2*(ndof+nb1+nb2)) = M2inv*reshape(A2tilde,[ndof,2*nb2]);
        end
        if nb1+nb2+1<=nb
            D(1:ndof,2*(ndof+nb1+nb2)+1:2*(ndof+nb1+nb2+nb3)) = M2inv*reshape(A3tilde,[ndof,2*nb3]);
        end
        if nb1+nb2+nb3+1<=nb
            D(1:ndof,2*(ndof+nb1+nb2+nb3)+1:2*(ndof+nb1+nb2+nb3+nb4)) = M2inv*reshape(A4tilde,[ndof,2*nb4]);
        end
    end
    D(ndof+1:2*ndof,1:ndof) = eye(2);

    aux = eye(2);
    for j=1:nb1-1
        aux = cat(1,aux,eye(2));
    end
    D(2*ndof+1:2*(ndof+nb1),3:4) = aux;
    sindex = 2*ndof+1;
    findex = 2*(ndof+nb1);
    array = sindex:findex;
    beta = repelem(B(1:nb1),2);
    for j=1:length(array)
        D(array(j),array(j)) = -U(i).*beta(j);
    end
    if nb1+1<=nb
        sindex = 2*(ndof+nb1)+1;
        findex = 2*(ndof+nb1+nb2);
        array = sindex:findex;
        D(sindex:findex,sindex-ndof:findex-ndof) = eye(2*nb2);
        beta = repelem(B(nb1+1:nb1+nb2),2);
        for j=1:length(array)
            D(array(j),array(j)) = -U(i).*beta(j);
        end
    end
    if nb1+nb2+1<=nb
        sindex = 2*(ndof+nb1+nb2)+1;
        findex = 2*(ndof+nb1+nb2+nb3);
        array = sindex:findex;
        D(sindex:findex,sindex-ndof:findex-ndof) = eye(2*nb3);
        beta = repelem(B(nb1+nb2+1:nb1+nb2+nb3),2);
        for j=1:length(array)
            D(array(j),array(j)) = -U(i).*beta(j);
        end
    end
    if nb1+nb2+nb3+1<=nb
        sindex = 2*(ndof+nb1+nb2+nb3)+1;
        findex = 2*(ndof+nb1+nb2+nb3+nb4);
        array = sindex:findex;
        D(sindex:findex,sindex-ndof:findex-ndof) = eye(2*nb4);
        beta = repelem(B(nb1+nb2+nb3+1:nb1+nb2+nb3+nb4),2);
        for j=1:length(array)
            D(array(j),array(j)) = -U(i).*beta(j);
        end
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

    % Velocidade Crítica de Flutter (Interpolação quadrática)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flutter_i_U==1
        eigval_ant = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
        eigval_pos = flutter_sol(U==U(flutter_i_U+2)).eigval(flutter_id);
        x = [real(eigval_flutter), real(eigval_ant), real(eigval_pos)];
        y = [U(flutter_i_U), U(flutter_i_U+1), U(flutter_i_U+2)];
        coefficients = polyfit(x, y, 2);
    else
        eigval_ant = flutter_sol(U==U(flutter_i_U-1)).eigval(flutter_id);
        eigval_pos = flutter_sol(U==U(flutter_i_U+1)).eigval(flutter_id);
        x = [real(eigval_ant), real(eigval_flutter), real(eigval_pos)];
        y = [U(flutter_i_U-1), U(flutter_i_U), U(flutter_i_U+1)];
        coefficients = polyfit(x, y, 2);
    end
    Uf_intquad = polyval(coefficients, 0);
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

        eigval_ant = flutter_sol(flutter_i_U-1).eigval(id);
        eigval_ft = flutter_sol(flutter_i_U).eigval(id);
        eigval_pos = flutter_sol(flutter_i_U+1).eigval(id-1);
        x = [imag(eigval_ant), imag(eigval_ft), imag(eigval_pos)];
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
[sizeS,dum] = size(S);
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
filename = sprintf('3_output/stability_analysis/optimized_%sForm/%s_poles/',...
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
    sprintf('2_pipeline/stability_analysis/optimized_%sForm/%s_poles/',...
    form,num2str(Np))];
movefile('*.dat',filename)