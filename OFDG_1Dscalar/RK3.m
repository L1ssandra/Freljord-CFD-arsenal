% RK3.m
t = 0;
T = 0;

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

while t < tend
    % calculate the damping term
    calculate_sigma
    
    dt = CFL/(max(abs(uh(:,1)))/hx + sigmaMAX/hx);
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    % Stage I
    du = Lh(uh);
    uh1 = uh + dt*du;
%     uh1 = TVD_Limiter(uh1);
    
    % Stage II
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
%     uh2 = TVD_Limiter(uh2);
    
    % Stage III
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
%     uh = TVD_Limiter(uh);
    
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
    
    fprintf('%d  %d  %d\n',t,max(abs(uh(:,1))),sigmaMAX)
     
end