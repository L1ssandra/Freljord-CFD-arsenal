% calculate_sigma.m
uhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1);
uhxR = zeros(Nx + 1,1);
uhxL = zeros(Nx + 1,1);
uhxxR = zeros(Nx + 1,1);
uhxxL = zeros(Nx + 1,1);
ujump = zeros(Nx + 1,1);
uxjump = zeros(Nx + 1,1);
uxxjump = zeros(Nx + 1,1);
uhb = [[0,0,0];uh;[0,0,0]];

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

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
    
    ujump(i) = uR - uL;
    uxjump(i) = uxR - uxL;
    uxxjump(i) = uxxR - uxxL;
end
sigmaMAX = 0;
for i = 1:Nx
    T1 = (1/9)*sqrt(ujump(i)^2 + ujump(i + 1)^2)/hx;
    T2 = (1/3)*sqrt(uxjump(i)^2 + uxjump(i + 1)^2)/hx;
    
    T3 = (4/135)*sqrt(ujump(i)^2 + ujump(i + 1)^2)/hx;
    T4 = (4/45)*sqrt(uxjump(i)^2 + uxjump(i + 1)^2)/hx;
    T5 = (2/27)*sqrt(uxxjump(i)^2 + uxxjump(i + 1)^2)/hx;
    
    sigmaMAX = max([T1,T2,T3,T4,T5,sigmaMAX]);
end