% main.m
clear;clc

global Nx dimPk NumGLP NumEq
Nx = 800;
k = 2;
dimPk = k + 1;
NumGLP = k + 1;

NumEq = 3;
CFL = 0.1;

get_GLP

init_data

get_basis

RK3

calculate_L2_Error

% draw_solution
draw_Leblanc
