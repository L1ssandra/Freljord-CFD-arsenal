% auto.m
clear;clc
global Nx Limit_type
Nx = 600;
Limit_type = 2;
main
u1 = uh;
Limit_type = 5;
main
u2 = uh;
Limit_type = 7;
main
u3 = uh;
Xc1 = Xc;

% calculate the "real" solution
Limit_type = 2;
Nx = 1500;
main
draw_solution1