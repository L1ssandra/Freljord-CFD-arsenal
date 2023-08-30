% init_data.m
global bcL bcR hx hx1 Xc lambda eta

%----------------
xa = 0;
xb = 2*pi;
u0 = @(x) sin(x); % Linear
% u0 = @(x) 0.5 + sin(x); % Burgers
eta = 0.005;
urealf = @(x,t) exp(-eta*t)*sin(x - t);
bcL = 1;
bcR = 1;
tend = 2;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

ureal = zeros(Nx,NumGLP);
for i = 1:Nx
    for j = 1:NumGLP
        ureal(i,j) = u0(Xc(i) + hx1*lambda(j));
    end
end

Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end