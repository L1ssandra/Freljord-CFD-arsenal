% main.m
%clear;clc

global Nx dimPk NumGLP %ujump uxjump uxxjump
Nx = 200;
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

drawaverage