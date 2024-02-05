% RK3.m
alpha = 1;
for i = 1:Nx
    for j = 1:NumGLP
        [SR,SL] = wavespeed(uh(i,j,:));
        alpha1 = max(abs(SR),abs(SL));
        if alpha1 > alpha
            alpha = alpha1;
        end
    end
end
dt = CFL*hx/alpha;
t = 0;

while t < tend
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    % Stage I
    du = Lh(uh);
    uh1 = uh + dt*du;
    uh1 = TVB_Limiter(uh1);
    
    % Stage II
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
    uh2 = TVB_Limiter(uh2);
    
    % Stage III
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
    uh = TVB_Limiter(uh);
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end