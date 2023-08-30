% drive the exact solution of Burgers' equation:
% u_t + u*u_x = 0, 0 < x < 2*pi, u(x,0) = a + b*sin(x),  ...(1)
% we can first solve u1_t + u1*u1_x = 0, u1(x,0) = sin(x),  ...(2)
% and let u(x,t) = a + b*u1(x - at,bt).

% clear;clc
% xa = -pi;
% xb = pi;
a = 0.5; b = 1; t = tend;
%Nx = 500;
Nx1 = Nx;% + 1;
hx = (xb - xa)/Nx;
hx1 = 0.5*hx;
X = xa + 0.5*hx:hx:xb - 0.5*hx;
X1 = zeros(1,Nx1);
for i = 1:Nx1
    X1((i - 1)*NumGLP + 1:i*NumGLP) = X(i) + hx1*lambda;
end
X = X1;
Nx1 = NumGLP*Nx;
u = zeros(1,Nx1);
epsilon = 1e-14;

u0 = @(x) sin(x);
du0 = @(x) cos(x);
method = 1;

tic
% solve (2) with x = x - at, t = bt, then we can get u1(x - at,bt)
X2 = X - a*t;
t2 = b*t;
% for each (x,t), we need to solve the equation of z:
% F(z) := z + u0(z)*t - x = 0
if method == 1
%% Bisect
% if 2k*pi < x < (2k + 1)*pi, then 
% F(2k*pi) = 0 + sin(0)*t - x < 0; F((2k + 1)*pi) = pi + sin(pi) - x > 0
% else, likewise, we have F(pi) < 0, F(2*pi) > 0.
    for i = 1:Nx1
        x = X2(i);
        F = @(z) z + u0(z)*t2 - x;
        delta = mod(x,2*pi);
        if mod(x,2*pi) < pi
            xa = x - delta;
            xb = xa + pi;
        else
            xa = x - delta + pi;
            xb = xa + pi;
        end
        xstar = bisect(F,xa,xb,epsilon);
        u(i) = u0(xstar);
    end
elseif method == 2
%% Newton
% we have F'(z) = 1 + u0'(z)*t
% the iteration function is G(z) = z - F(z)/F'(z)
% for each (x,t), if x < pi, we take x0 = 0
%                    x > pi, we take x0 = 2*pi
    for i = 1:Nx1
        x = X2(i);
        F = @(z) z + u0(z)*t2 - x;
        dF = @(z) 1 + du0(z)*t2;
        G = @(z) z - F(z)/dF(z);
        delta = mod(x,2*pi);
        if mod(x,2*pi) < pi
            x0 = x - delta;
        else
            x0 = x - delta + 2*pi;
        end
        xstar = fpi(G,x0,epsilon,100);
        u(i) = u0(xstar);
    end
end
time = toc;

% change to the solution of (1)
u = a + b*u;

figure(1)
hold on
plot(X,u,'b-','linewidth',1.3)
% plot(pi + a*t,a,'rx')
% axis([X(1),X(end),min(u) - 0.05,max(u) + 0.05]);