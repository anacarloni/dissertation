clearvars
close all
clc
addpath('figures')
addpath('functions')

% Identificação dos ROMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idERA
idBPOD
idOKIDERA
idVolterra

% Resposta ao Impulso a partir de Unit Sample (US) and Discrete Step (DS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
impRespUSandDS

% Resposta dos ROMs aos Sinais WF4 e DS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arbitraryInputWF4
arbitraryInputDS

% Resposta ao Impulso dos ROMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idImpulseResp
impulseResp

% Análise de Estabilidade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROMid