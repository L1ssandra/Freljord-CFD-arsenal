% Spectral-DG 2D
% u_t + f(u)_x + g(u)_y = 0, u(x,y,0) = sin(cos(x + y))
% x - direction: DG
% y - direction: Spectral
%clear;clc

global Nx dimPk NumGLP L
Nx = 60; % DG
k = 2; % DG polynomial degree
NumGLP = 5;
dimPk = k + 1;
L = 4; % Spectral
CFL = 0.15;

get_GLP

init_data

get_basis

L2_Pro

RK3

%flash2D

calculate_L2_error
