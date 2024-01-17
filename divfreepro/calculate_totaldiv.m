% calculate_totaldiv.m
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

dBxG = zeros(Nx,Ny1,NumGLP);
for i = 1:Nx
    for j = 1:Ny
        for d = 1:dimPk
            for i1 = 1:NumGLP
                dBxG(i,j,i1) = dBxG(i,j,i1) + uh(i,j,d,6)*phixG(i1,d);
            end
        end
    end
end

divG = dBxG + dByG;
totaldiv = 0; Qdiv = divG(:,:,3);
for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            totaldiv = totaldiv + hx1*weight(i1)*abs(divG(i,j,i1));
        end
    end
end