% calculate_L2error_order.m
clear;clc
global Nx
AccNum = 6;
Table = zeros(AccNum,4);
Nx0 = 20;

for nn = 1:AccNum
    Nx = Nx0*2^(nn - 1);
    main
    Table(nn,1) = L1_Error;
    Table(nn,3) = L2_Error;
    Table(nn,5) = L8_Error;
    if nn > 1
        Table(nn,2) = log(Table(nn - 1,1)/Table(nn,1))/log(2);
        Table(nn,4) = log(Table(nn - 1,3)/Table(nn,3))/log(2);
        Table(nn,6) = log(Table(nn - 1,5)/Table(nn,5))/log(2);
    end
end