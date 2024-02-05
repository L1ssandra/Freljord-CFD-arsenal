% calculate_L2_Error.m

uE = abs(uh - ureal);
L2_Error = 0;
L1_Error = 0;

% L2_Error
for i = 1:Nx
    for i1 = 1:NumGLP
        L2_Error = L2_Error + hx1*weight(i1)*uE(i,i1,1)^2;
    end
end
L2_Error = sqrt(L2_Error);

% L1_Error
for i = 1:Nx
    for i1 = 1:NumGLP
        L1_Error = L1_Error + hx1*weight(i1)*abs(uE(i,i1,1));
    end
end