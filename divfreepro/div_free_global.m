% div_free_global.m
c = zeros(Nx,Ny,dimPk1);

uhG7 = zeros(Nx,Ny1,NumGLP);
dByG = zeros(Nx,Ny1,NumGLP);
for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for d = 1:dimPk1
                uhG7(i,j,i1) = uhG7(i,j,i1) + uh(i,j,d,7)*phiG(i1,d);
            end
        end
    end
end

% calculate the term d(By)/dy
for i = 1:Nx
    for i1 = 1:NumGLP
        
        gusin = zeros(1,L);
        gucos = zeros(1,L);
        guysin = zeros(1,L);
        guycos = zeros(1,L);
        for j = 1:Ny
            for d = 1:L
                gusin(d) = gusin(d) + sinkY(d,j)*uhG7(i,j,i1);
                gucos(d) = gucos(d) + coskY(d,j)*uhG7(i,j,i1);
            end
        end
        gusin = gusin/L;
        gucos = gucos/L;
        
        for d = 1:L
            guysin(d) = -d*gucos(d);
            guycos(d) = d*gusin(d);
        end
        
        for j = 1:Ny
            for d = 1:L
                dByG(i,j,i1) = dByG(i,j,i1) + guysin(d)*sinkY(d,j) + guycos(d)*coskY(d,j);
            end
        end
        
    end
end

% project d(By)/dy to x-axis
for i = 1:Nx
    for j = 1:Ny
        for d = 1:dimPk1
            for i1 = 1:NumGLP
                c(i,j,d) = c(i,j,d) + 0.5*weight(i1)*dByG(i,j,i1)*phiG(i1,d)/mm(d);
            end
        end
    end
end
c0M = -c(:,:,1); c1M = -c(:,:,2); c2M = -c(:,:,3);

Bxedge = zeros(Nx + 1,Ny);
for j = 1:Ny
    Bxedge(1,j) = sum(uh(:,j,1,6))/Nx;
    for i = 1:Nx
        Bxedge(i + 1,j) = Bxedge(i,j) + hx*c0M(i,j);
    end
end

uh1 = zeros(Nx,Ny1,dimPk,NumEq);
% modify Bx
for i = 1:Nx
    for j = 1:Ny
        e0 = c0M(i,j); e1 = c1M(i,j); e2 = c2M(i,j);
        d3 = hx/6*e2; d2 = hx/4*e1; d1 = 3*d3/5 + hx/2*(e0 - e2/3);
        d0 = 0.5*(Bxedge(i,j) + Bxedge(i + 1,j)) - 2/3*d2;
        uh1(i,j,1,6) = d0;
        uh1(i,j,2,6) = d1;
        uh1(i,j,3,6) = d2;
        uh1(i,j,4,6) = d3;
    end
end
uh1bar = zeros(1,Ny);
uhbar = zeros(1,Ny);
for j = 1:Ny
uhbar(j) = sum(uh(:,j,1,6))/Nx;
uh1bar(j) = sum(uh1(:,j,1,6))/Nx;
uh(:,j,1,6) = uh1(:,j,1,6) + uhbar(j) - uh1bar(j);
end
uh(:,:,2:dimPk,6) = uh1(:,:,2:dimPk,6);

