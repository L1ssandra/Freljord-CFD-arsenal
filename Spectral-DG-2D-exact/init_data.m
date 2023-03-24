% init_data.m
global bcL bcR hx hx1 hy f g Ny Ny1 R Xc Y lambda

%----------------
xa = 0;
xb = 2*pi;
ya = 0;
yb = 2*pi;
%u0 = @(x,y) sin(cos(x + y));
%u0 = @(x,y) exp(-((x - pi).^2 + (y - pi).^2));
u0 = @(x,y) sin(x)*sin(y);
f = @(u) u; g = @(u) u;
R = @(x,y,t) sin(x - t).*cos(y);
bcL = 1;
bcR = 1;
tend = 1;
%----------------

Ny = 2*L; Ny1 = Ny + 1;
hx = (xb - xa)/Nx; hy = (yb - ya)/Ny;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;
Y = ya:hy:yb;

ureal = zeros(Nx,Ny1,NumGLP);
for i = 1:Nx
    for j = 1:Ny1
        for i1 = 1:NumGLP
            ureal(i,j,i1) = u0(Xc(i) + hx1*lambda(i1),Y(j));
        end
    end
end