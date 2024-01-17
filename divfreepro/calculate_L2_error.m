uhG = zeros(Nx,Ny1,NumGLP,NumEq);

for i = 1:Nx
    for j = 1:Ny1
        for d = 1:dimPk
            for i1 = 1:NumGLP
                for n = 1:NumEq
                    uhG(i,j,i1,n) = uhG(i,j,i1,n) + uh(i,j,d,n)*phiG(i1,d);
                end
            end
        end
    end
end

uE = abs(uhG - ureal);
L2_Error = zeros(1,8);

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for n = 1:NumEq
                L2_Error(n) = L2_Error(n) + 0.5*weight(i1)*uE(i,j,i1,n)^2;
            end
        end
    end
end
L2_Error = sqrt(L2_Error/(Nx*Ny));

% L8_Error = max(max(max(abs(uE))));