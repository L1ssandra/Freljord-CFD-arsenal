function uh = TVB_Limiter(uh)

global Nx bcL bcR Mvec NumGLP NumEq gamma

uhb = zeros(Nx + 2,NumGLP,NumEq);
uhb(1,:,:) = zeros(1,NumGLP,NumEq);
uhb(end,:,:) = zeros(1,NumGLP,NumEq);
uhb(2:end - 1,:,:) = uh;

uibar = zeros(Nx + 2,NumEq);
uhmod = zeros(Nx,NumGLP,NumEq);

deltaUR = zeros(NumEq,1); deltaUL = zeros(NumEq,1);
deltaURM = zeros(NumEq,1); deltaULM = zeros(NumEq,1);
deltaURM1 = zeros(NumEq,1); deltaULM1 = zeros(NumEq,1);
thetai = zeros(NumEq,1);

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

% calculate the cell-average
for i = 1:Nx + 2
    for n = 1:NumEq
        uibar(i,n) = 0.5*uhb(i,:,n)*Mvec;
    end
end

for i = 1:Nx

    for n = 1:NumEq
        deltaUR(n) = uh(i,NumGLP,n) - uibar(i + 1,n);
        deltaUL(n) = uibar(i + 1,n) - uh(i,1,n);
        deltaURM(n) = uibar(i + 2,n) - uibar(i + 1,n);
        deltaULM(n) = uibar(i + 1,n) - uibar(i,n);
    end
    
    v = uibar(i + 1,2)/uibar(i + 1,1);
    p = pressure(uibar(i + 1,1),uibar(i + 1,2),uibar(i + 1,3));
    c = sqrt(gamma*p/uibar(i + 1,1));
    H = (uibar(i + 1,3) + p)/uibar(i + 1,1);
    
    R1 = [1;v - c;H - v*c];
    R2 = [1;v;0.5*v^2];
    R3 = [1;v + c;H + v*c];
    
    R = [R1,R2,R3];%+ 1e-12*eye(3);
    
    L = inv(R);
    
    deltaURM1 = R*minmod(L*deltaUR,L*deltaURM,L*deltaULM);
    deltaULM1 = R*minmod(L*deltaUL,L*deltaURM,L*deltaULM);
    
    for n = 1:NumEq
        uhmod(i,NumGLP,n) = uibar(i + 1,n) + deltaURM1(n);
        uhmod(i,1,n) = uibar(i + 1,n) - deltaULM1(n);
    end
    
    for n = 1:NumEq
        thetai(n) = (deltaURM1(n) - deltaULM1(n))/(deltaUR(n) - deltaUL(n) + 1e-16);
    end
    
    for d = 2:NumGLP - 1
        for n = 1:NumEq
            uhmod(i,d,n) = uibar(i + 1,n) + thetai(n)*(uh(i,d,n) - uibar(i + 1,n));
        end
    end
end

uh = uhmod;

end

function a1 = minmod(a,b,c)

global hx M NumEq

a1 = zeros(NumEq,1);

for i = 1:NumEq
    if abs(a(i)) < M*hx^2
        a1(i) = a(i);
    else
        if sign(a(i)) == sign(b(i)) && sign(a(i)) == sign(c(i))
            a1(i) = sign(a(i))*min(abs([a(i),b(i),c(i)]));
        else
            a1(i) = 0;
        end
    end
end

end