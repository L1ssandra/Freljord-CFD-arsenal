% IIF3.m
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
    epsilon = 1e-2;
    
    if IIFindex < 3
        dt = epsilon*hx^2;
    else
        dt = CFL*hx/alpha;
    end
    
    if t + dt >= tend
        dt = tend - t;
        t = tend;
    else
        t = t + dt;
    end
    
    du = Lh(uh);
    
    duIIF(end,:) = du;
    dtIIF(end) = dt;
    
%     Cmat = construct_Cmat(uh);
    Cmat = construct_Cmat_4th(uh);
    
    dt1 = dtIIF(2); dt2 = dtIIF(1);
    
    beta3 = 1 + 1/(dt1*(dt1 + dt2))*( dt^2/3 + dt/2*(2*dt1 + dt2) ); 
    beta2 = -1/(dt1*dt2)*(dt^2/3 + dt/2*(dt1 + dt2)); 
    beta1 = 1/(dt2*(dt1 + dt2))*(dt^2/3 + dt*dt1/2);
    
    beta2_2nd = 1/dt1*(dt/2 + dt1); 
    beta1_2nd = -dt/(2*dt1);
    
    if IIFindex == 1 % IIF1
        uh = (expm(dt*Cmat)*uh')' + dt*(expm(dt*Cmat)*duIIF(3,:)')';
%         uh = exp(-t)*sin(Xc - t);
    elseif IIFindex == 2 % IIF2
        uh = (expm(dt*Cmat)*uh')' + dt*(beta2_2nd*expm(dt*Cmat)*duIIF(3,:)' + beta1_2nd*expm((dt + dt1)*Cmat)*duIIF(2,:)')';
%         uh = exp(-t)*sin(Xc - t);
    else % IIF3
        uh = (expm(dt*Cmat)*uh')' + dt*(beta3*expm(dt*Cmat)*duIIF(3,:)' + beta2*expm((dt + dt1)*Cmat)*duIIF(2,:)' + beta1*expm((dt + dt1 + dt2)*Cmat)*duIIF(1,:)')';
    end
    
    IIFindex = min([IIFindex + 1,3]);
    
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