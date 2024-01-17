% get_basis.m
global phiG phixG mm phiGR phiGL sinkY coskY

% DG
phiG = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
    phiG(i,4) = lambda(i)^3 - 3/5*lambda(i);
    
    phixG(i,1) = 0;
    phixG(i,2) = 1/hx1;
    phixG(i,3) = 2*lambda(i)/hx1;
    phixG(i,4) = (3*lambda(i)^2 - 3/5)/hx1;
end

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;
phiGR(4) = 2/5;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;
phiGL(4) = -2/5;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;
mm(4) = 4/175;

% Spectral
sinkY = zeros(L,Ny1);
coskY = zeros(L,Ny1);
for d = 1:L
    sinkY(d,:) = sin(d*Y);
    coskY(d,:) = cos(d*Y);
end