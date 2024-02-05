function du = Lh(uh)

global Nx NumGLP bcL bcR NumEq
global Minv B S hx1
uhb = zeros(Nx + 2,NumGLP,NumEq);
uhb(1,:,:) = zeros(1,NumGLP,NumEq);
uhb(end,:,:) = zeros(1,NumGLP,NumEq);
uhb(2:end - 1,:,:) = uh;
uhR = zeros(Nx + 1,NumEq);
uhL = zeros(Nx + 1,NumEq);
fhat = zeros(Nx + 1,NumEq);
du = zeros(Nx,NumGLP,NumEq);
fvec = zeros(Nx,NumGLP,NumEq);

% set_bc
if bcL == 1
    uhb(1,:,:) = uh(end,:,:);
elseif bcL == 2
    uhb(1,:,:) = uh(2,:,:);
end

if bcR == 1
    uhb(end,:,:) = uh(1,:,:);
elseif bcR == 2
    uhb(end,:,:) = uh(end - 1,:,:);
end

% Step 1: calculate the f in cell
for i = 1:Nx
    for j = 1:NumGLP
        fvec(i,j,1) = f1(uh(i,j,1),uh(i,j,2),uh(i,j,3));
        fvec(i,j,2) = f2(uh(i,j,1),uh(i,j,2),uh(i,j,3));
        fvec(i,j,3) = f3(uh(i,j,1),uh(i,j,2),uh(i,j,3));
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for n = 1:NumEq
        uhR(i,n) = uhb(i,NumGLP,n);
        uhL(i,n) = uhb(i + 1,1,n);
    end
end

for i = 1:Nx + 1
    uR = uhL(i,:);
    uL = uhR(i,:);
    fR = [f1(uR(1),uR(2),uR(3)),f2(uR(1),uR(2),uR(3)),f3(uR(1),uR(2),uR(3))];
    fL = [f1(uL(1),uL(2),uL(3)),f2(uL(1),uL(2),uL(3)),f3(uL(1),uL(2),uL(3))];
    [SR_R,SL_R] = wavespeed(uR);
    [SR_L,SL_L] = wavespeed(uL);
    SR = max(SR_R,SR_L); SL = min(SL_R,SL_L);
    
    fhat(i,:) = HLL_Flux(uR,uL,fR,fL,SR,SL);
end

% Step 3: calculate du
for i = 1:Nx
    for n = 1:NumEq
        fhatveclocal = [fhat(i,n);zeros(NumGLP - 2,1);fhat(i + 1,n)];
        fveclocal = fvec(i,:,n)';
        % weak form
        % du(i,:) = Minv*(S'*fveclocal - B*fhatveclocal)/hx1;
        % strong form
        du(i,:,n) = Minv*(B*(fveclocal - fhatveclocal) - S*fveclocal)/hx1;
    end
end

end