uhG = zeros(Nx,Ny1,NumGLP);

for i = 1:Nx
    for j = 1:Ny1
        for d = 1:dimPk
            for i1 = 1:NumGLP
                uhG(i,j,i1) = uhG(i,j,i1) + uh(i,j,d)*phiG(i1,d);
            end
        end
    end
end

uE = abs(uhG - ureal);
L2_Error = 0;

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            L2_Error = L2_Error + 0.5*weight(i1)*uE(i,j,i1)^2;
        end
    end
end
L2_Error = sqrt(L2_Error/(Nx*Ny));

L8_Error = max(max(max(abs(uE))));