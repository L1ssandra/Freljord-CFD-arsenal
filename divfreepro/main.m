% Spectral-DG 2D: MHD equation
% x - direction: DG
% y - direction: Spectral, the solution should be periodic along y-direction

clear;clc

global Nx dimPk NumGLP L NumEq dimPk1
Nx = 100; % DG
k = 2; % DG polynomial degree
NumEq = 8; % 2.5-dimensional MHD
NumGLP = 5;
dimPk = k + 2;
dimPk1 = k + 1;
L = 8; % Spectral
CFL = 0.15;

get_GLP

init_data

get_basis

L2_Pro

% div_free_local
div_free_global

% RK3

% flash2D

calculate_L2_error

% fprintf('%d  %d\n',L2_Error(6),L2_Error(7))

calculate_totaldiv

C0test

fprintf('||f_x + g_y|| = %d\n',totaldiv)
fprintf('||[uh]|| = %d\n',jump)
