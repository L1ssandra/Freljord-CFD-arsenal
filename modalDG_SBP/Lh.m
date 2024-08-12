function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx hx1
global Pmat Lmat DNmat Nfmat
uhb = [[0,0,0];uh;[0,0,0]];
uhG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);
fint = zeros(Nx,NumGLP + 2);

% uh(:,2:end) = 0;

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% Step 1: calculate the Integral in cell
for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for d = 1:dimPk
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);
    end
end

for i = 1:Nx
    for i1 = 1:NumGLP
        fint(i,i1) = f(uhG(i,i1));
    end
    fint(i,NumGLP + 1) = f(uhL(i + 1));
    fint(i,NumGLP + 2) = f(uhR(i));
end

for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    alpha = max(abs(df(uR)),abs(df(uL)));
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
end

for i = 1:Nx
    fvec = fint(i,:)';
%     size([Pmat,Lmat])
%     size([DNmat])
%     size(fint)
    du(i,:) = -[Pmat,Lmat]*DNmat*fvec - Lmat*Nfmat*(  [fhat(i + 1);fhat(i)] - fvec(NumGLP + 1:NumGLP + 2)  );
end

du = du/hx1;

% du(:,2:end) = 0;

end
