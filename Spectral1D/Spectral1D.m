% Spectral1D.m
% spectral method solving ut + f(u)x = 0

global f sinkX coskX L Nx fuh Nx1 uh uc
Nplot = 80;
xa = 0;
xb = 2*pi;
L = 5;
Nx = 2*L;
Nx1 = 2*L + 1;
Nplot1 = Nplot + 1;
hx = (xb - xa)/Nx;
hplot = (xb - xa)/Nplot;
X = xa:hx:xb;
Xplot = xa:hplot:xb;
CFL = 0.1;
tend = 2*pi;
t = 0;
T = 0;

RKorder = 3;

u0 = @(x) sin(x);
f = @(x) x;

uh = zeros(1,Nx1);
fuh = zeros(1,Nx1);

sinkX = zeros(L,Nx1);
coskX = zeros(L,Nx1);

sinkXplot = zeros(L,Nplot1);
coskXplot = zeros(L,Nplot1);

for d = 1:L
    sinkX(d,:) = sin(d*X);
    coskX(d,:) = cos(d*X);
    sinkXplot(d,:) = sin(d*Xplot);
    coskXplot(d,:) = cos(d*Xplot);
end

uc = 0;
usin = zeros(1,L);
ucos = zeros(1,L);

fuc = 0;
fusin = zeros(1,L);
fucos = zeros(1,L);

duc = 0;
dusin = zeros(1,L);
ducos = zeros(1,L);

for i = 1:Nx1
    uh(i) = u0(X(i));
    fuh(i) = f(uh(i));
end

% Fourier transform
for i = 1:Nx
    uc = uc + uh(i);
    for d = 1:L
        usin(d) = usin(d) + sinkX(d,i)*uh(i);
        ucos(d) = ucos(d) + coskX(d,i)*uh(i);
    end
    
end
uc = uc/Nx;
usin = usin/L;
ucos = ucos/L;

uplot = uc*ones(1,Nplot1);

for i = 1:Nplot1
    for d = 1:L
        uplot(i) = uplot(i) + sinkXplot(d,i)*usin(d) + coskXplot(d,i)*ucos(d);
    end
end

% plot(Xplot,uplot,'r.',Xplot,u0(Xplot),'b-');
% axis([Xplot(1),Xplot(end),min(uplot) - 0.1,max(uplot) + 0.1]);
uflash = uplot;

while t < tend
    
    alpha = 1;
    dt = CFL*hx/alpha;
    
    if t + dt < tend
        t = t + dt;
    else
        dt = tend - t;
        t = tend;
    end
    
    if RKorder == 1 % Euler-Forward
        [dusin,ducos] = Lh(usin,ucos);
        usin = usin + dt*dusin;
        ucos = ucos + dt*ducos;
    elseif RKorder == 3 % RK3
        % Stage I
        [dusin,ducos] = Lh(usin,ucos);
        usin1 = usin + dt*dusin;
        ucos1 = ucos + dt*ducos;
        % Stage II
        [dusin,ducos] = Lh(usin1,ucos1);
        usin2 = (3/4)*usin + (1/4)*usin1 + (1/4)*dt*dusin;
        ucos2 = (3/4)*ucos + (1/4)*ucos1 + (1/4)*dt*ducos;
        % Stage III
        [dusin,ducos] = Lh(usin2,ucos2);
        usin = (1/3)*usin + (2/3)*usin2 + (2/3)*dt*dusin;
        ucos = (1/3)*ucos + (2/3)*ucos2 + (2/3)*dt*ducos;
    end
        
    uplot = uc*ones(1,Nplot1);

    for i = 1:Nplot1
        for d = 1:L
            uplot(i) = uplot(i) + sinkXplot(d,i)*usin(d) + coskXplot(d,i)*ucos(d);
        end
    end
    T = [T;t];
    uflash(end + 1,:) = uplot;
    
end

L2error = sqrt(sum((uplot - u0(Xplot)).^2)/Nplot1);
L8error = max(uplot - u0(Xplot));
