% get_basis.m
global phiG phiGLL phixG mm phiGR phiGL phixxGL phixxGR
global phixGL phixGR
phiG = zeros(NumGLP,dimPk);
phiGLL = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phixxGL = zeros(1,dimPk);
phixxGR = zeros(1,dimPk);
phixGL = zeros(1,dimPk);
phixGR = zeros(1,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
    
    phiGLL(i,1) = 1;
    phiGLL(i,2) = lambdaL(i);
    phiGLL(i,3) = lambdaL(i)^2 - 1/3;
    
    phixG(i,1) = 0;
    phixG(i,2) = 1/hx1;
    phixG(i,3) = 2*lambda(i)/hx1;
    
end

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;

phixGR(1) = 0;
phixGR(2) = 1/hx1;
phixGR(3) = 2/hx1;

phixGL(1) = 0;
phixGL(2) = 1/hx1;
phixGL(3) = -2/hx1;

phixxGR(1) = 0;
phixxGR(2) = 0;
phixxGR(3) = 2/hx1^2;

phixxGL = phixxGR;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;