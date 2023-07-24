function [J,Fc,kFc,A] = fiterror_2ndForm(Bini,H,M,k,kk)

% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

% Multiplicidade dos Polos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[B,p,~] = multiplicity(Bini);

% Parâmetros da Multiplicidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb = size(B,2);                 % número total de polos iniciais
nb1 = p(1);                     % índices de polos não repetidos
nb2 = p(1)+p(2);                % índices de polos com multiplicidade 2
nb3 = p(1)+p(2)+p(3);           % índices de polos com multiplicidade 3
N2 = nb-nb1;
N3 = 2*nb-nb1-nb2;
N4 = 3*nb-nb1-nb2-nb3;
N = size(k,1);
nn = size(kk,1);
Np = p(1)+2*p(2)+3*p(3)+4*p(4);
n = Np+3;

% Strings de Frequencia (Eversman e Tewari - 2a e 3a formas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc = zeros(N,3+nb);
kFc = zeros(nn,3+nb);
for i = 1:N
    Fc(i,1:3+nb) = [1 1i*k(i) -(k(i)^2) 1./(1i*k(i)+B(1:nb))];
    if nb1+1<=nb
        Fc(i,1:3+nb+N2) = [Fc(i,1:3+nb) 1./(1i*k(i)+B(nb1+1:nb)).^2];
    end
    if nb2+1<=nb
        Fc(i,1:3+nb+N3) = [Fc(i,1:3+nb+N2) 1./(1i*k(i)+B(nb2+1:nb)).^3];
    end
    if nb3+1<=nb
        Fc(i,1:3+nb+N4) = [Fc(i,1:3+nb+N3) 1./(1i*k(i)+B(nb3+1:nb)).^4];
    end
end

for i = 1:nn
    kFc(i,1:3+nb) = [1 1i*kk(i) -(kk(i)^2) 1./(1i*kk(i)+B(1:nb))];
    if nb1+1<=nb
        kFc(i,1:3+nb+N2) = [kFc(i,1:3+nb) 1./(1i*kk(i)+B(nb1+1:nb)).^2];
    end
    if nb2+1<=nb
        kFc(i,1:3+nb+N3) = [kFc(i,1:3+nb+N2) 1./(1i*kk(i)+B(nb2+1:nb)).^3];
    end
    if nb3+1<=nb
        kFc(i,1:3+nb+N4) = [kFc(i,1:3+nb+N3) 1./(1i*kk(i)+B(nb3+1:nb)).^4];
    end
end

% Inicialização de Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = struct('Clh',zeros(n,1),'Cla',zeros(n,1),...
    'Cmh',zeros(n,1),'Cma',zeros(n,1));
epsilon = struct('Clh',zeros(1),'Cla',zeros(1),...
    'Cmh',zeros(1),'Cma',zeros(1));

% Eversman e Tewari (1991)
% Pseudo-Inversa de Penrose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AUX = Fc.'*conj(Fc) + Fc'*Fc;
AUX1 = AUX\eye(size(AUX));

A.Clh = AUX1*(Fc'*H.Clh + Fc.'*conj(H.Clh));
A.Cla = AUX1*(Fc'*H.Cla + Fc.'*conj(H.Cla));
A.Cmh = AUX1*(Fc'*H.Cmh + Fc.'*conj(H.Cmh));
A.Cma = AUX1*(Fc'*H.Cma + Fc.'*conj(H.Cma));

% Erro das Aproximações
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon.Clh = (abs(Fc*A.Clh - H.Clh).^2)./M.Clh;
epsilon.Cla = (abs(Fc*A.Cla - H.Cla).^2)./M.Cla;
epsilon.Cmh = (abs(Fc*A.Cmh - H.Cmh).^2)./M.Cmh;
epsilon.Cma = (abs(Fc*A.Cma - H.Cma).^2)./M.Cma;

% Função Objetivo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = sum([epsilon.Clh epsilon.Cla epsilon.Cmh epsilon.Cma],'all');

% Restrição de Desigualdade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nb
    if B(i)<=0
        J = 1e30;
    end
end

end