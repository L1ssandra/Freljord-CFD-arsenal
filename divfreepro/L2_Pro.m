% L2_Pro.m

uh = zeros(Nx,Ny1,dimPk,NumEq);

for i = 1:Nx
    for j = 1:Ny1
        for d = 1:dimPk1
            for i1 = 1:NumGLP
                for n = 1:NumEq
                    uh(i,j,d,n) = uh(i,j,d,n) + 0.5*weight(i1)*ureal(i,j,i1,n)*phiG(i1,d);
                end
            end
        end
    end
end

for d = 1:dimPk
    uh(:,:,d,:) = uh(:,:,d,:)/mm(d);
end