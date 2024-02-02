% main.m
% clear;clc

global Nx dimPk NumGLP
% Nx = 100;
k = 2;
dimPk = k + 1;
NumGLP = k + 1;
CFL = 0.2;

get_GLP

init_data

get_basis

RK3

calculate_L2_Error

draw_solution
