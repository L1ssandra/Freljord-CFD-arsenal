% IIF1.m

t = 0;
T = 0;

frameMAX = 200;
uhflash = zeros(1,Nx);

uhflash(1,:) = uh;
t1 = tend/frameMAX;
ii1 = 1;

while t < tend
    
    alpha = 1;
    
    dt = CFL*hx/alpha;
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    du = Lh(uh);
    
%     Cmat = construct_Cmat(uh);
    Cmat = construct_Cmat_4th(uh);
    
    uh = (expm(dt*Cmat)*uh')' + dt*du;
    
    if t >= ii1*t1
        uhflash(end + 1,:) = uh;
        T = [T;t];
        ii1 = ii1 + 1;
    end
    
    fprintf('%d  %d\n',t,max(abs(uh(:,1))))
     
end