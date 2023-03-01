% L2_Pro.m

uh = zeros(Nx,Ny1,dimPk);

for i = 1:Nx
    for j = 1:Ny1
        for d = 1:dimPk
            for i1 = 1:NumGLP
                uh(i,j,d) = uh(i,j,d) + 0.5*weight(i1)*ureal(i,j,i1)*phiG(i1,d);
            end
        end
    end
end

for d = 1:dimPk
    uh(:,:,d) = uh(:,:,d)/mm(d);
end