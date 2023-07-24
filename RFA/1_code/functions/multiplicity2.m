function [B,mult,p,index] = multiplicity2(Bi)

% Comparação dos Polos Bi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np = size(Bi,2);
Baux = zeros(Np,2);
threshold = .01;
for i=1:Np-1
    cont = 1;
    for j=i+1:Np
        if abs(Bi(i)-Bi(j)) <= threshold
            cont = cont+1;
%             Bi(j) = Bi(i);
        end
    end
    Baux(i,:) = [Bi(i) cont];
end
Baux(Np,:) = [Bi(Np) 1];

Baux2 = zeros(Np,2);
for i=1:Np-1
    cont = 1;
    for j=i+1:Np
        if abs(Bi(i)-Bi(j)) <= threshold
            cont = cont+1;
            Bi(j) = Bi(i);
        end
    end
    Baux2(i,:) = [Bi(i) cont];
end
Baux2(Np,:) = [Bi(Np) 1];

% % Polos Reduzidos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,~,ic] = uniquetol(Baux(:,1),threshold);
% a_counts = accumarray(ic,1);
% value = [C,a_counts];
% % B = value(:,1)';
% % idout = value(:,2)';
% 
% % Polos Removidos
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~,rmvid] = unique(ic,'first');
% for i=1:Np
%     if isempty(find(rmvid==i))==1
%         rmv(i) = i;
%     end
% end
% rmv = find(rmv~=0);
% 
% [value(:,2),id] = sort(value(:,2));
% value(:,1) = value(id,1);
% B = value(:,1)';
% 
% [val,id] = unique(Baux(:,1),'stable');
% B = val';

% Multiplicidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,id] = unique(Baux2(:,1),'stable');
mult = Baux2(id,2)';
p = zeros(1,4);
p(1) = length(find(mult==1));
p(2) = length(find(mult==2));
p(3) = length(find(mult==3));
p(4) = length(find(mult==4));

% Ordena os Polos com Multiplicidade Crescente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = .01;
[C,IA,IC] = uniquetol(Bi,threshold);
[id,ord] = sort(IC,'ascend');
B = Bi(ord);

% [tblB,index] =sortrows(Baux);

Bflip = flip(Baux,2);
[tblBflip,index] =sortrows(Bflip);
Bord = flip(tblBflip,2);
B = Bord(:,1);
mult = tblBflip(:,1);

end