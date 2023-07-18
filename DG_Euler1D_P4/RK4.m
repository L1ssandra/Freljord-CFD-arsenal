% SSPRK4.m

CFL = 0.6;
dt = CFL*hx;
t = 0;

uh = TVD_Limiter(uh);

while t < tend
    
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
        uh1 = TVD_Limiter(uh1);
    end
    
    uh2 = 0.04*uh2 + 0.36*uh1;
    uh1 = 15*uh2 - 5*uh1;
    
    % Stage II
    for i = 6:9
        du = Lh(uh1);
        uh1 = uh1 + (1/6)*dt*du;
        uh1 = TVD_Limiter(uh1);
    end
    
    % Stage III
    du = Lh(uh1);
    uh = uh2 + 0.6*uh1 + 0.1*dt*du;
    uh = TVD_Limiter(uh);
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end