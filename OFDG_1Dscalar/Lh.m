function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx
global phixGR phixGL phixxGR phixxGL s
% global ujump uxjump uxxjump
uhb = [[0,0,0];uh;[0,0,0]];
uhG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
uhxR = zeros(Nx + 1,1);
uhxL = zeros(Nx + 1,1);
uhxxR = zeros(Nx + 1,1);
uhxxL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
ujump = zeros(Nx + 1,1);
uxjump = zeros(Nx + 1,1);
uxxjump = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);

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

for i = 1:Nx
    for d = 2:dimPk
        for i1 = 1:NumGLP
            du(i,d) = du(i,d) + 0.5*weight(i1)*f(uhG(i,i1))*phixG(i1,d);
        end
    end
end

% Step 2: calculate the flux at edge
for i = 1:Nx + 1
    for d = 1:dimPk
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);
        uhxR(i) = uhxR(i) + uhb(i,d)*phixGR(d);
        uhxL(i) = uhxL(i) + uhb(i + 1,d)*phixGL(d);
        uhxxR(i) = uhxxR(i) + uhb(i,d)*phixxGR(d);
        uhxxL(i) = uhxxL(i) + uhb(i + 1,d)*phixxGL(d);
    end
end

for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    uxR = uhxL(i);
    uxL = uhxR(i);
    uxxR = uhxxL(i);
    uxxL = uhxxR(i);
    
    alpha = max(abs(uR),abs(uL));
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
    
    ujump(i) = uR - uL;
    uxjump(i) = uxR - uxL;
    uxxjump(i) = uxxR - uxxL;
end

% Step 3: calculate [f*phi] on each cell
for i = 1:Nx
    for d = 1:dimPk
        du(i,d) = du(i,d) - (1/hx)*(phiGR(d)*fhat(i + 1) - phiGL(d)*fhat(i));
    end
end

% Step 4: add the damping term
for i = 1:Nx
    du(i,2) = du(i,2) - s*(1/9)*uh(i,2)*sqrt(ujump(i)^2 + ujump(i + 1)^2)/hx;
    du(i,2) = du(i,2) - s*(1/3)*uh(i,2)*sqrt(uxjump(i)^2 + uxjump(i + 1)^2)/hx;
    
    du(i,3) = du(i,3) - s*(4/135)*uh(i,3)*sqrt(ujump(i)^2 + ujump(i + 1)^2)/hx;
    du(i,3) = du(i,3) - s*(4/45)*uh(i,3)*sqrt(uxjump(i)^2 + uxjump(i + 1)^2)/hx;
    du(i,3) = du(i,3) - s*(2/27)*uh(i,3)*sqrt(uxxjump(i)^2 + uxxjump(i + 1)^2)/hx;
end

for d = 1:dimPk
    du(:,d) = du(:,d)/mm(d);
end

