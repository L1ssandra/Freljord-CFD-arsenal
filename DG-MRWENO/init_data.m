% init_data.m
global bcL bcR hx hx1 M s Ck gammaC gammaN

%----------------
xa = 0;
xb = 2*pi;
% u0 = @(x) 0.5 + sin(x);
bcL = 1;
bcR = 1;
tend = 20*pi;
M = 1;
s = 3;
Ck = 10;
gammaC = 10/11;
gammaN = 1/11;
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