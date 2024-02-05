% init_data.m
global bcL bcR hx hx1 Xc M gamma

%----------------
gamma = 1.4;

xa = -10;
xb = 10;

% rho0 = @(x) 1 + 0.2*sin(x);
% u0 = @(x) 1 + 0.*x;
% p0 = @(x) 1 + 0.*x;

% rho0 = @(x) 1 + 0.*x;
% u0 = @(x) 0.*x;
% p0 = @(x) 1 + 0.*x;

rho0 = @(x) tube_rho0(x);
u0 = @(x) tube_u0(x);
p0 = @(x) tube_p0(x);

bcL = 2;
bcR = 2;
tend = 0.0001;

M = 1;
%----------------

U1 = @(x) rho0(x);
U2 = @(x) rho0(x).*u0(x);
U3 = @(x) p0(x)./(gamma - 1) + 0.5*rho0(x).*u0(x).*u0(x);
% U3 = @(x) 1e-12 + 0.*x; % Blast

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
uh = ureal;
% uh(Nx/2,:,3) = 3200000/hx; % Blast