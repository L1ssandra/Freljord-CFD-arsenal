% C0test.m
BxR = zeros(Nx,Ny1);
BxL = zeros(Nx,Ny1);
for i = 1:Nx
    for j = 1:Ny
        for d = 1:dimPk
            BxR(i,j) = BxR(i,j) + uh(i,j,d,6)*phiGR(d);
            BxL(i,j) = BxL(i,j) + uh(i,j,d,6)*phiGL(d);
        end
    end
end

jump = 0;
for i = 1:Nx - 1
    for j = 1:Ny
        jump = jump + abs(BxL(i + 1,j) - BxR(i,j))/Nx;
    end
end
        