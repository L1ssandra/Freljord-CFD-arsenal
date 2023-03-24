% RK3.m
global tRK
dt = CFL/(1/hx + 1/hy);
t = 0; T = 0;
fprintf('%d  %d\n',t,max(max(abs(uh(:,:,1)))))
uflash = uh(:,:,1);

while t < tend
    
    tRK = t;
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    % Stage I
    du = Lh(uh);
    uh1 = uh + dt*du;
    
    % Stage II
    tRK = tRK + dt;
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
    
    % Stage III
    tRK = tRK - 0.5*dt;
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
    
    fprintf('%d  %d\n',t,max(max(abs(uh(:,:,1)))))
    uflash(:,:,end + 1) = uh(:,:,1);
    T = [T;t];
     
end