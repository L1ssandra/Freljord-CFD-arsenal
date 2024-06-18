% init_data.m
global bcL bcR hx hx1 Xc
%----------------
xa = 0;
xb = 2*pi;
u0 = @(x) sin(x);
bcL = 1;
bcR = 1;
tend = 1;
%----------------

hx = (xb - xa)/Nx;
hx1 = 0.5*hx;

Xc = xa + hx1:hx:xb - hx1;

uh = zeros(1,Nx);
ureal = zeros(1,Nx);
for i = 1:Nx
    uh(i) = u0(Xc(i));
    ureal(i) = exp(-tend)*u0(Xc(i) - tend);
%     ureal(i) = u0(Xc(i) - tend);
end
