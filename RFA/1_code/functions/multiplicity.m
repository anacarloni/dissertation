function [B,p,idout] = multiplicity(Bi)

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
            Bi(j) = Bi(i);
        end
    end
    Baux(i,:) = [Bi(i) cont];
end
Baux(Np,:) = [Bi(Np) 1];

% Polos Reduzidos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C,~,ic] = unique(Baux(:,1));
a_counts = accumarray(ic,1);
value = [C,a_counts];
B = value(:,1)';
idout = value(:,2)';

% Polos Removidos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [~,rmvid] = unique(ic,'first');
% for i=1:Np
%     if isempty(find(rmvid==i))==1
%         rmv(i) = i;
%     end
% end
% rmv = find(rmv~=0);

% [value(:,2),id] = sort(value(:,2));
% value(:,1) = value(id,1);
% B = value(:,1)';

% [val,id] = unique(Baux(:,1),'stable');
% B = val';

% Multiplicidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,id] = unique(Baux(:,1),'stable');
mult = Baux(id,2)';
p = zeros(1,4);
p(1) = length(find(mult==1));
p(2) = length(find(mult==2));
p(3) = length(find(mult==3));
p(4) = length(find(mult==4));

end