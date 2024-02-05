function uh = TVB_Limiter(uh)

global Nx dimPk bcL bcR Mvec NumGLP

uhb = [zeros(1,NumGLP);uh;zeros(1,NumGLP)];
uibar = zeros(1,Nx + 2);
uhmod = zeros(Nx,dimPk);
uhmod(:,1) = uh(:,1);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% calculate the cell-average 
for i = 1:Nx + 2
    uibar(i) = 0.5*uhb(i,:)*Mvec;
end

for i = 1:Nx
    deltaUR = uh(i,NumGLP) - uibar(i + 1);
    deltaUL = uibar(i + 1) - uh(i,1);
    deltaURM = uibar(i + 2) - uibar(i + 1);
    deltaULM = uibar(i + 1) - uibar(i);
    
    deltaURM1 = minmod(deltaUR,deltaURM,deltaULM);
    deltaULM1 = minmod(deltaUL,deltaURM,deltaULM);
    
    uhmod(i,NumGLP) = uibar(i + 1) + deltaURM1;
    uhmod(i,1) = uibar(i + 1) - deltaULM1;
    
    thetai = (deltaURM1 - deltaULM1)/(deltaUR - deltaUL + 1e-16);
    
    for d = 2:NumGLP - 1
        uhmod(i,d) = uibar(i + 1) + thetai*(uh(i,d) - uibar(i + 1));
    end
end

uh = uhmod;

end

function a1 = minmod(a,b,c)

global hx M

if abs(a) < M*hx^2
    a1 = a;
else
    if sign(a) == sign(b) && sign(a) == sign(c)
        a1 = sign(a)*min(abs([a,b,c]));
    else
        a1 = 0;
    end
end

end