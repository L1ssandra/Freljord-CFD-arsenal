% main.m
% clear;clc

global Nx dimPk NumGLP
% Nx = 100;
k = 2;
NumGLP = 5;
dimPk = k + 1;
CFL = 0.2;

get_GLP

init_data

get_basis

L2_Pro

RK3

calculate_L2_Error

% plot(Xc,uh(:,1),'-b','linewidth',1); axis([0,2*pi,-0.1,1.1])
% plot(Xc,uh(:,1),'-b','linewidth',1); axis([0,2*pi,-1.1,1.1])
% flash
