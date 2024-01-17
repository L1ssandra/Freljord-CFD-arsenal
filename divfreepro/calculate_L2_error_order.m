% calculate_L2_error_order.m

Nk = 6;
Table = zeros(Nk,2);
Nx0 = 20;

for nn = 1:Nk
    Nx = Nx0*2^(nn - 1);
    main
    Table(nn,1) = L2_Error(6);
    Table(nn,3) = L2_Error(7);
%     Table(n,3) = L8_Error;
    if nn > 1
        Table(nn,2) = log(Table(nn - 1,1)/Table(nn,1))/log(2);
        Table(nn,4) = log(Table(nn - 1,3)/Table(nn,3))/log(2);
%         Table(n,4) = log(Table(n - 1,3)/Table(n,3))/log(2);
    end
end