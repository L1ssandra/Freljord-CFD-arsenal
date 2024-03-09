% init_data.m
global bcL bcR hx hx1 M gamma

%----------------
xa = -0.5;
xb = 0.5;

gamma = 1.4;

rho = @(x) Sod_rho0(x);%1 + 0.2*sin(x);%
u = @(x) Sod_u0(x);%1 + 0.*x;%
p = @(x) Sod_p0(x);%1 + 0.*x;%

U1 = @(x) rho(x); % rho
U2 = @(x) rho(x).*u(x); % rho u
U3 = @(x) p(x)./(gamma - 1) + 0.5*rho(x).*u(x).^2; % E
% U3 = @(x) 1e-12 + 0.*x; %E

M = 1;
bcL = 2;
bcR = 2;
tend = 0.2;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

ureal = zeros(Nx,NumGLP,NumEq);
for i = 1:Nx
    for j = 1:NumGLP
        ureal(i,j,1) = U1(Xc(i) + hx1*lambda(j));
        ureal(i,j,2) = U2(Xc(i) + hx1*lambda(j));
        ureal(i,j,3) = U3(Xc(i) + hx1*lambda(j));
    end
end