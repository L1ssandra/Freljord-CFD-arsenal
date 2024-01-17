% init_data.m
global bcL bcR hx hx1 hy f g Ny Ny1 gamma gamma1

%----------------
gamma = 5/3; gamma1 = gamma - 1;
xa = 0;
xb = 2*pi;
ya = 0;
yb = 2*pi;
rho0 = @(x,y) 1;
u0 = @(x,y) 1;
v0 = @(x,y) 1;
w0 = @(x,y) 1;
p0 = @(x,y) 1;
Bx0 = @(x,y) 4 + sin(x).*sin(y);
By0 = @(x,y) 4 + cos(x).*cos(y);
Bz0 = @(x,y) 0;

rhou0 = @(x,y) rho0(x,y).*u0(x,y);
rhov0 = @(x,y) rho0(x,y).*v0(x,y);
rhow0 = @(x,y) rho0(x,y).*w0(x,y);
E0 = @(x,y) p0(x,y)/gamma1 + 0.5*rho0(x,y).*(u0(x,y).^2 + v0(x,y).^2 + w0(x,y).^2)  + 0.5*(Bx0(x,y).^2 + By0(x,y).^2 + Bz0(x,y).^2);

U1 = @(x,y) rho0(x,y);
U2 = @(x,y) rhou0(x,y);
U3 = @(x,y) rhov0(x,y);
U4 = @(x,y) rhow0(x,y);
U5 = @(x,y) E0(x,y);
U6 = @(x,y) Bx0(x,y);
U7 = @(x,y) By0(x,y);
U8 = @(x,y) Bz0(x,y);

f = @(u) u; g = @(u) u;
bcL = 1;
bcR = 1;
tend = 0;
%----------------

Ny = 2*L; Ny1 = Ny + 1;
hx = (xb - xa)/Nx; hy = (yb - ya)/Ny;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;
Y = ya:hy:yb;

ureal = zeros(Nx,Ny1,NumGLP,NumEq);
for i = 1:Nx
    for j = 1:Ny1
        for i1 = 1:NumGLP
            ureal(i,j,i1,1) = U1(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,2) = U2(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,3) = U3(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,4) = U4(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,5) = U5(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,6) = U6(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,7) = U7(Xc(i) + hx1*lambda(i1),Y(j));
            ureal(i,j,i1,8) = U8(Xc(i) + hx1*lambda(i1),Y(j));
        end
    end
end