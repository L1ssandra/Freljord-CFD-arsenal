function du = Lh(uh)

global Nx NumGLP dimPk bcL bcR
global Minv B S hx1
uhb = [zeros(1,NumGLP);uh;zeros(1,NumGLP)];
uhG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);
fvec = zeros(Nx,NumGLP); fveclocal = zeros(NumGLP,1);
fhatveclocal = zeros(NumGLP,1);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% Step 1: calculate the f in cell
for i = 1:Nx
    for j = 1:NumGLP
        fvec(i,j) = f(uh(i,j));
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    uhR(i) = uhb(i,NumGLP);
    uhL(i) = uhb(i + 1,1);
end

for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    alpha = max(abs(uR),abs(uL));
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
end

% Step 3: calculate du
for i = 1:Nx
    fhatveclocal = [fhat(i);zeros(NumGLP - 2,1);fhat(i + 1)];
    fveclocal = fvec(i,:)';
    % weak form
    % du(i,:) = Minv*(S'*fveclocal - B*fhatveclocal)/hx1;
    % strong form
    du(i,:) = Minv*(B*(fveclocal - fhatveclocal) - S*fveclocal)/hx1;
end

end