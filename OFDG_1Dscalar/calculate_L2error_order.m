% calculate_L2error_order.m

Table = zeros(7,4);
Nx0 = 20;

for n = 1:7
    Nx = Nx0*2^(n - 1);
    main
    Table(n,1) = L2_Error;
    Table(n,3) = L8_Error;
    if n > 1
        Table(n,2) = log(Table(n - 1,1)/Table(n,1))/log(2);
        Table(n,4) = log(Table(n - 1,3)/Table(n,3))/log(2);
    end
end