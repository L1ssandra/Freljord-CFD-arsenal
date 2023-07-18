% main.m
%clear;clc

global Nx dimPk NumGLP NumEq flux_type
Nx = 200;
k = 4;
NumGLP = 5;
dimPk = k + 1;
NumEq = 3;
CFL = 0.1;
flux_type = 3;
Limit_type = 4;

get_GLP

init_data

get_basis
 
L2_Pro
 
RK3

calculate_L2_Error

draw_solution
