% RK4.m
frameMAX = 200;
uhflash = zeros(1,Nx*NumGLP);

uhG = zeros(Nx,NumGLP);

for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end
uhG = reshape(uhG',Nx*NumGLP,1)';

uhflash(1,:) = uhG;
t1 = tend/frameMAX;
i1 = 1;

T = 0;
t = 0;

while t < tend
    
    alphaM = 1;
    for i = 1:Nx
        if abs(dH(uh(i,1))) > alphaM
            alphaM = abs(dH(uh(i,1)));
        end
    end
    
    dt = CFL*hx/alphaM;
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    uh1 = uh;
    uh2 = uh;
    
    % Stage I
    for i = 1:5
        du = Lh(uh1);
        uh1 = uh1 + (1/6)*dt*du;
    end
    
    uh2 = 0.04*uh2 + 0.36*uh1;
    uh1 = 15*uh2 - 5*uh1;
    
    % Stage II
    for i = 6:9
        du = Lh(uh1);
        uh1 = uh1 + (1/6)*dt*du;
    end
    
    % Stage III
    du = Lh(uh1);
    uh = uh2 + 0.6*uh1 + 0.1*dt*du;
    
    if t >= i1*t1
        uhG = zeros(Nx,NumGLP);
        for i = 1:Nx
            for d = 1:dimPk
                uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
            end
        end
        uhG = reshape(uhG',Nx*NumGLP,1)';
        uhflash(end + 1,:) = uhG;
        T = [T;t];
        i1 = i1 + 1;
    end
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end