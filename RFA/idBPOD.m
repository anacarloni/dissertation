clearvars
close all
clc

%% Identificar ROM a partir do ERA usando BPOD

% Rank, r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=2:10

    % ROM identificado usando ERA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([fileparts(pwd),filesep,'2_pipeline/idERA/sysERAlowAmp_r' num2str(r) '.mat']);
    sysERA = sysERAlowAmp;

    % Identificar matrizes do ROM usando BPOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rBPOD =  r;
    [~,~,xFull] = impulse(sysERA,0:1:(rBPOD*5)+1);
    sysAdjoint = ss(sysERA.A',sysERA.C',sysERA.B',sysERA.D',-1);
    [yAdjoint,~,xAdjoint] = impulse(sysAdjoint,0:1:(rBPOD*5)+1);
    HankelOC = [];  % compute Hankel matrix H=OC
    for i=2:size(xAdjoint,1) % we start at 2 to avoid incorporating the D matrix
        Hrow = [];
        for j=2:size(xFull,1)
            Ystar = permute(squeeze(xAdjoint(i,:,:)),[2 1]);
            MarkovParameter = Ystar*squeeze(xFull(j,:,:));
            Hrow = [Hrow MarkovParameter];
        end
        HankelOC = [HankelOC; Hrow];
    end
    [U,Sig,V] = svd(HankelOC);
    Xdata = [];
    Ydata = [];
    for i=2:size(xFull,1)  % we start at 2 to avoid incorporating the D matrix
        Xdata = [Xdata squeeze(xFull(i,:,:))];
        Ydata = [Ydata squeeze(xAdjoint(i,:,:))];
    end
    Phi = Xdata*V*Sig^(-1/2); % modes
    Psi = Ydata*U*Sig^(-1/2);
    Ar = Psi(:,1:rBPOD)'*sysERA.a*Phi(:,1:rBPOD);
    Br = Psi(:,1:rBPOD)'*sysERA.b;
    Cr = sysERA.c*Phi(:,1:rBPOD);
    Dr = sysERA.d;
    sysBPOD = ss(Ar,Br,Cr,Dr,-1);

    % Salvar o ROM identificado com BPOD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename = sprintf('2_pipeline/idBPOD/sysBPOD_r%d.mat',r);
    save([fileparts(pwd),filesep,filename],'sysBPOD');
end