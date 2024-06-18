% RK3.m
% global t dt
t = 0;
T = 0;

frameMAX = 200;
uhflash = zeros(1,Nx);

uhflash(1,:) = uh;
t1 = tend/frameMAX;
ii1 = 1;

duIIF = zeros(3,Nx);
dtIIF = zeros(1,3);
IIFindex = 1;

% beta3 = 23/12; beta2 = -4/3; beta1 = 5/12;
% beta2_2nd = 3/2; beta1_2nd = -1/2;

while t < tend
    
    alpha = 1;
    
    dt = CFL*hx^2;
    
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
    du = Lh(uh1);
    uh2 = (3/4)*uh + (1/4)*uh1 + (1/4)*dt*du;
    
    % Stage III
    du = Lh(uh2);
    uh = (1/3)*uh + (2/3)*uh2 + (2/3)*dt*du;
    
    if t >= ii1*t1
        uhflash(end + 1,:) = uh;
        T = [T;t];
        ii1 = ii1 + 1;
    end
    
    fprintf('%d  %d\n',t,max(abs(uh)))
    
end