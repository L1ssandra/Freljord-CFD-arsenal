% SSPRK4.m

CFL = 0.6;
t = 0;

while t < tend
    
    alpha = 1;
    for i = 1:Nx
        [alpha1,~] = wavespeed(uh(i,1,:));
        if alpha1 > alpha
            alpha = alpha1;
        end
    end
    dt = CFL*hx/alpha;
    
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
        if Limit_type == 2
            uh1 = TVD_Limiter_P2(uh1);
        elseif Limit_type == 3
            uh1 = TVD_Limiter_P3(uh1);
        elseif Limit_type == 5
            uh1 = SWENO_Limiter(uh1);
        elseif Limit_type == 6
            uh1 = SWENOZ_Limiter(uh1);
        elseif Limit_type == 7
            uh1 = MRWENO_Limiter(uh1);
        end
        uh1 = pp_Limiter(uh1);
    end
    
    uh2 = 0.04*uh2 + 0.36*uh1;
    uh1 = 15*uh2 - 5*uh1;
    
    % Stage II
    for i = 6:9
        du = Lh(uh1);
        uh1 = uh1 + (1/6)*dt*du;
        if Limit_type == 2
            uh1 = TVD_Limiter_P2(uh1);
        elseif Limit_type == 3
            uh1 = TVD_Limiter_P3(uh1);
        elseif Limit_type == 5
            uh1 = SWENO_Limiter(uh1);
        elseif Limit_type == 7
            uh1 = MRWENO_Limiter(uh1);
        end
        uh1 = pp_Limiter(uh1);
    end
    
    % Stage III
    du = Lh(uh1);
    uh = uh2 + 0.6*uh1 + 0.1*dt*du;
    if Limit_type == 2
        uh = TVD_Limiter_P2(uh);
    elseif Limit_type == 3
        uh = TVD_Limiter_P3(uh);
    elseif Limit_type == 5
        uh = SWENO_Limiter(uh);
    elseif Limit_type == 6
        uh = SWENOZ_Limiter(uh);
    elseif Limit_type == 7
        uh = MRWENO_Limiter(uh);
    end
    uh = pp_Limiter(uh);
    
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end