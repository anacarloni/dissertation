function [J,Fc,kFc,A] = fiterror_1stForm(B,H,M,k,kk)

% Inicialização
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

% Parâmetros da Multiplicidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb = size(B,2);
N = size(k,1);
nn = size(kk,1);
n = nb+3;

% Strings de Frequencia (Eversman e Tewari - 1a forma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fc = zeros(N,3+nb);
kFc = zeros(nn,3+nb);
for i = 1:N
    Fc(i,1:3+nb) = [1 1i*k(i) -(k(i)^2) 1./(1i*k(i)+B(1:nb))];
end

for i = 1:nn
    kFc(i,1:3+nb) = [1 1i*kk(i) -(kk(i)^2) 1./(1i*kk(i)+B(1:nb))];
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
    if B(i)<0
        J = 1e30;
    end
end

end