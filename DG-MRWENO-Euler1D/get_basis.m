% get_basis.m
global phiG phiGLL phixG mm phiGR phiGL
phiG = zeros(NumGLP,dimPk);
phiGLL = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
    phiG(i,4) = lambda(i)^3 - 3/5*lambda(i);
%     phiG(i,5) = lambda(i)^4 - 6/7*lambda(i)^2 + 3/35;
    
    phiGLL(i,1) = 1;
    phiGLL(i,2) = lambdaL(i);
    phiGLL(i,3) = lambdaL(i)^2 - 1/3;
    phiGLL(i,4) = lambdaL(i)^3 - 3/5*lambdaL(i);
%     phiGLL(i,5) = lambdaL(i)^4 - 6/7*lambdaL(i)^2 + 3/35;
    
    phixG(i,1) = 0;
    phixG(i,2) = 1/hx1;
    phixG(i,3) = 2*lambda(i)/hx1;
    phixG(i,4) = (3*lambda(i)^2 - 3/5)/hx1;
%     phixG(i,5) = (4*lambda(i)^3 - 12/7*lambda(i))/hx1;
end

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;
phiGR(4) = 2/5;
% phiGR(5) = 8/35;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;
phiGL(4) = -2/5;
% phiGL(5) = 8/35;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;
mm(4) = 4/175;
% mm(5) = 64/11025;