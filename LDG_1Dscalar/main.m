% main.m
% clear;clc

global Nx dimPk NumGLP
% Nx = 100;
k = 3;
NumGLP = 5;
dimPk = k + 1;
CFL = 1/(k^4);

get_GLP

init_data

get_basis

L2_Pro

% RK3
RK4

calculate_L2_Error

% drawaverage
draw_solution