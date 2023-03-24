function du = Lh(uh)
global Nx Ny Ny1 NumGLP dimPk mm phiG phixG phiGR phiGL
global bcL bcR weight hx f g R sinkY coskY L tRK Xc Y lambda hx1
uhb = zeros(Nx + 2,Ny1,dimPk);
uhb(2:end - 1,:,:) = uh;
uhG = zeros(Nx,Ny1,NumGLP);
fuG = zeros(Nx,Ny1,NumGLP);
guyG = zeros(Nx,Ny1,NumGLP);
uhR = zeros(Nx + 1,Ny1);
uhL = zeros(Nx + 1,Ny1);
fhat = zeros(Nx + 1,Ny1);
du = zeros(Nx,Ny1,dimPk);
RHS = zeros(Nx,Ny1,NumGLP);

% set_bc
if bcL == 1
    uhb(1,:,:) = uh(end,:,:); 
end

if bcR == 1
    uhb(end,:,:) = uh(1,:,:);
end

% Step 1: calculate the Integral in cell
for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for d = 1:dimPk
                uhG(i,j,i1) = uhG(i,j,i1) + uh(i,j,d)*phiG(i1,d);
            end
        end
    end
end

for i = 1:Nx
    for j = 1:Ny
        for d = 2:dimPk
            for i1 = 1:NumGLP
                du(i,j,d) = du(i,j,d) + 0.5*weight(i1)*f(uhG(i,j,i1))*phixG(i1,d);
            end
        end
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for j = 1:Ny
        for d = 1:dimPk
            uhR(i,j) = uhR(i,j) + uhb(i,j,d)*phiGR(d);
            uhL(i,j) = uhL(i,j) + uhb(i + 1,j,d)*phiGL(d);
        end
    end
end

for i = 1:Nx + 1
    for j = 1:Ny
        uR = uhL(i,j);
        uL = uhR(i,j);
        alpha = 1;
        fhat(i,j) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
    end
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx
    for j = 1:Ny
        for d = 1:dimPk
            du(i,j,d) = du(i,j,d) - (1/hx)*(phiGR(d)*fhat(i + 1,j) - phiGL(d)*fhat(i,j));
        end
    end
end

% Step 4: Fourier transform
for i = 1:Nx
    for i1 = 1:NumGLP
        gusin = zeros(1,L);
        gucos = zeros(1,L);
        guysin = zeros(1,L);
        guycos = zeros(1,L);
        for j = 1:Ny
            for d = 1:L
                gusin(d) = gusin(d) + sinkY(d,j)*g(uhG(i,j,i1));
                gucos(d) = gucos(d) + coskY(d,j)*g(uhG(i,j,i1));
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
                guyG(i,j,i1) = guyG(i,j,i1) + guysin(d)*sinkY(d,j) + guycos(d)*coskY(d,j);
            end
        end
        
        % source term
        for j = 1:Ny
            RHS(i,j,i1) = R(Xc(i) + hx1*lambda(i1),Y(j),tRK);
        end
        
    end
end

RHS = RHS - guyG;

for i = 1:Nx
    for j = 1:Ny
        for d = 1:dimPk
            for i1 = 1:NumGLP
                du(i,j,d) = du(i,j,d) + 0.5*weight(i1)*RHS(i,j,i1)*phiG(i1,d);
            end
        end
    end
end

for d = 1:dimPk
    du(:,:,d) = du(:,:,d)/mm(d);
end

du(:,end,:) = du(:,1,:);

end