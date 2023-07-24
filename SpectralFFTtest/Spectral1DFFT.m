%Spectral1DFFT.m
% spectral method solving ut + f(u)x = 0 using Fast Fourier Transform

global f sinkX coskX L Nx fuh Nx1 uh uc L1
% Nplot = 80;
xa = 0;
xb = 2*pi;
% L = 8;
L1 = L - 1;
Nx = 2*L;
Nx1 = 2*L + 1;
Nplot1 = Nplot + 1;
hx = (xb - xa)/Nx;
hplot = (xb - xa)/Nplot;
X = xa:hx:xb - hx;
Xplot = xa:hplot:xb;
CFL = 0.1;
tend = 20*pi;
t = 0;
T = 0;

RKorder = 3;

u0 = @(x) sin(cos(x));
f = @(x) x;

uh = zeros(1,Nx);
fuh = zeros(1,Nx);

sinkX = zeros(L1,Nx);
coskX = zeros(L1,Nx);

sinkXplot = zeros(L1,Nplot1);
coskXplot = zeros(L1,Nplot1);

for d = 1:L1
    sinkX(d,:) = sin(d*X);
    coskX(d,:) = cos(d*X);
    sinkXplot(d,:) = sin(d*Xplot);
    coskXplot(d,:) = cos(d*Xplot);
end

uc = 0;

fuc = 0;
fusin = zeros(1,L1);
fucos = zeros(1,L1);

duc = 0;
dusin = zeros(1,L1);
ducos = zeros(1,L1);

for i = 1:Nx
    uh(i) = u0(X(i));
    fuh(i) = f(uh(i));
end

% Fourier transform
for i = 1:Nx
    uc = uc + uh(i);
end
uc = uc/Nx;
[ucos,usin] = pointtocoeff(uh);

uplot = uc*ones(1,Nplot1);

for i = 1:Nplot1
    for d = 1:L1
        uplot(i) = uplot(i) + sinkXplot(d,i)*usin(d) + coskXplot(d,i)*ucos(d);
    end
end

% plot(Xplot,uplot,'r.',Xplot,u0(Xplot),'b-');
% axis([Xplot(1),Xplot(end),min(uplot) - 0.1,max(uplot) + 0.1]);
uflash = uplot;

usin1 = zeros(1,L1);
usin2 = zeros(1,L1);
ucos1 = zeros(1,L1);
ucos2 = zeros(1,L1);

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
        [dusin,ducos] = LhFFT(usin,ucos);
        usin = usin + dt*dusin;
        ucos = ucos + dt*ducos;
    elseif RKorder == 3 % RK3
        % Stage I
        [dusin,ducos] = LhFFT(usin,ucos);
        usin1 = usin + dt*dusin;
        ucos1 = ucos + dt*ducos;
        % Stage II
        [dusin,ducos] = LhFFT(usin1,ucos1);
        usin2 = (3/4)*usin + (1/4)*usin1 + (1/4)*dt*dusin;
        ucos2 = (3/4)*ucos + (1/4)*ucos1 + (1/4)*dt*ducos;
        % Stage III
        [dusin,ducos] = LhFFT(usin2,ucos2);
        usin = (1/3)*usin + (2/3)*usin2 + (2/3)*dt*dusin;
        ucos = (1/3)*ucos + (2/3)*ucos2 + (2/3)*dt*ducos;
    end
        
    uplot = uc*ones(1,Nplot1);

    for i = 1:Nplot1
        for d = 1:L1
            uplot(i) = uplot(i) + sinkXplot(d,i)*usin(d) + coskXplot(d,i)*ucos(d);
        end
    end
    T = [T;t];
    uflash(end + 1,:) = uplot;
    
end

L2error = sqrt(sum((uplot - u0(Xplot)).^2)/Nplot1);
L8error = max(uplot - u0(Xplot));