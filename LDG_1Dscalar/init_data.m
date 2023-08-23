% init_data.m
global bcL bcR hx hx1

%----------------
xa = 0;
xb = 2*pi;
u0 = @(x) sin(x);
urealf = @(x,t) exp(-t)*sin(x - t);
bcL = 1;
bcR = 1;
tend = 0.01;
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