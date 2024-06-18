% IIF2.m

t = 0;
T = 0;

frameMAX = 200;
uhflash = zeros(1,Nx);

uhflash(1,:) = uh;
t1 = tend/frameMAX;
ii1 = 1;

duIIF = zeros(2,Nx);
dtIIF = zeros(1,2);
IIFindex = 1;

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
    
    duIIF(end,:) = du;
    dtIIF(end) = dt;
    
    dt1 = dtIIF(1);
    
    beta2 = 1/dt1*(dt/2 + dt1); 
    beta1 = -dt/(2*dt1);
    
%     Cmat = construct_Cmat(uh);
    Cmat = construct_Cmat_4th(uh);
    
    if IIFindex == 1 % IIF1
        uh = (expm(dt*Cmat)*uh')' + dt*(expm(dt*Cmat)*duIIF(2,:)')';
    else % IIF2
        uh = (expm(dt*Cmat)*uh')' + dt*(beta2*expm(dt*Cmat)*duIIF(2,:)' + beta1*expm((dt + dt1)*Cmat)*duIIF(1,:)')';
    end
    
    IIFindex = min([IIFindex + 1,2]);
    
    if t >= ii1*t1
        uhflash(end + 1,:) = uh;
        T = [T;t];
        ii1 = ii1 + 1;
    end
    
    fprintf('%d  %d\n',t,max(abs(uh)))
    
    % update the solution
    duIIF(1:end - 1,:) = duIIF(2:end,:);
    dtIIF(1:end - 1) = dtIIF(2:end);
     
end