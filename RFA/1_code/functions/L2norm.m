% Aquisição de Dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename1 = sprintf('0_data/fit_in_a_c00_case%s.dat',caso);
filename2 = sprintf('0_data/fit_in_h_c00_case%s.dat',caso);
fileID1 = fopen([fileparts(pwd),filesep,filename1]);
fileID2 = fopen([fileparts(pwd),filesep,filename2]);
C1 = textscan(fileID1,'%f %f %f %f %f');
C2 = textscan(fileID2,'%f %f %f %f %f');
fclose(fileID1);
fclose(fileID2);
[~,recla_c00,imcla_c00,recma_c00,imcma_c00] = deal(C1{:,1:5});
[k_c00,reclh_c00,imclh_c00,recmh_c00,imcmh_c00] = deal(C2{:,1:5});
N_c00 = size(k_c00,1);
kk_c00 = linspace(0,max(k_c00),1001)';

% Coeficientes Aerodinâmicos do Case 00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clhC00 = complex(reclh_c00,imclh_c00);
claC00 = complex(recla_c00,imcla_c00);
cmhC00 = complex(recmh_c00,imcmh_c00);
cmaC00 = complex(recma_c00,imcma_c00);

% Coeficientes Aerodinâmicos da Walsh Function (WF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clhWF = complex(reclh,imclh);
claWF = complex(recla,imcla);
cmhWF = complex(recmh,imcmh);
cmaWF = complex(recma,imcma);

% Coeficientes Aerodinâmicos do RFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clhRFA = CLhfit;
claRFA = CLafit;
cmhRFA = CMhfit;
cmaRFA = CMafit;

% Norma L2 entre Case 00 e WF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clhL2_WF = norm(clhWF-clhC00);
claL2_WF = norm(claWF-claC00);
cmhL2_WF = norm(cmhWF-cmhC00);
cmaL2_WF = norm(cmaWF-cmaC00);

% Norma L2 entre Case 00 e RFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clhL2_RFA = norm(clhRFA-clhC00);
claL2_RFA = norm(claRFA-claC00);
cmhL2_RFA = norm(cmhRFA-cmhC00);
cmaL2_RFA = norm(cmaRFA-cmaC00);
