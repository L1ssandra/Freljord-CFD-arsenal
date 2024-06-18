% calculate_L2_Error.m
uE = abs(uh - ureal);
L1_Error = 0;
L2_Error = 0;
L8_Error = max(max(uE));

for i = 1:Nx
    L1_Error = L1_Error + uE(i);
    L2_Error = L2_Error + uE(i)^2;
end
L1_Error = L1_Error/Nx;
L2_Error = sqrt(L2_Error/Nx);