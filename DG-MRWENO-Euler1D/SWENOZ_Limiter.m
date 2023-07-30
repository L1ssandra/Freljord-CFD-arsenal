function uh = SWENOZ_Limiter(uh)
% S = Simple
global dimPk Nx bcL bcR phiGL phiGR hx Ck gammaC gammaN NumEq gamma
uhb = zeros(Nx + 2,dimPk,NumEq);
uhb(2:end - 1,:,:) = uh;
uhR = zeros(Nx + 1,NumEq); uhL = zeros(Nx + 1,NumEq);
ujump = zeros(Nx + 1,NumEq);

uhmod = uh;
epsilon = 1e-10;

% set_bc
if bcL == 1
    uhb(1,:,:) = uh(end,:,:);
elseif bcL == 2
    uhb(1,:,:) = uh(2,:,:);
end

if bcR == 1
    uhb(end,:,:) = uh(1,:,:);
elseif bcR == 2
    uhb(end,:,:) = uh(end - 1,:,:);
end

% KXRCF Detector: find the trouble cells
for i = 1:Nx + 1
    for n = 1:NumEq
        for d = 1:dimPk
            uhR(i,n) = uhR(i,n) + uhb(i,d,n)*phiGR(d);
            uhL(i,n) = uhL(i,n) + uhb(i + 1,d,n)*phiGL(d);
        end
    end
end
% calculate the jump
for i = 1:Nx + 1
    ujump(i,:) = uhL(i,:) - uhR(i,:);
end
% calculate the indicator
is_trouble_cell = zeros(1,Nx);
Ind = zeros(1,Nx);
for i = 1:Nx
    for n = 1:NumEq
        Ind = (abs(ujump(i + 1,n)) + abs(ujump(i,n)))/(hx^2*max([abs(uhL(i,n)),abs(uhR(i,n))]) + 1e-16);
        if Ind > Ck
            is_trouble_cell(i) = 1;
        end
    end
end

for i = 1:Nx
    if is_trouble_cell(i) == 1
        polx = zeros(NumEq,1);
        polxx = zeros(NumEq,1);
        polxxx = zeros(NumEq,1);
        
        polxL = zeros(NumEq,1);
        polxxL = zeros(NumEq,1);
        polxxxL = zeros(NumEq,1);
        
        polxR = zeros(NumEq,1);
        polxxR = zeros(NumEq,1);
        polxxxR = zeros(NumEq,1);
        for n = 1:NumEq
            polx(n,1) = uhb(i + 1,2,n);
            polxx(n,1) = uhb(i + 1,3,n);
            polxxx(n,1) = uhb(i + 1,4,n);
            
            polxL(n,1) = uhb(i,2,n);
            polxxL(n,1) = uhb(i,3,n);
            polxxxL(n,1) = uhb(i,4,n);
            
            polxR(n,1) = uhb(i + 2,2,n);
            polxxR(n,1) = uhb(i + 2,3,n);
            polxxxR(n,1) = uhb(i + 2,4,n);
        end
        
        v = uh(i,1,2)/uh(i,1,1);
        p = pressure(uh(i,1,1),uh(i,1,2),uh(i,1,3));
        c = sqrt(gamma*p/uh(i,1,1));
        H = (uh(i,1,3) + p)/uh(i,1,1);

        R1 = [1;v - c;H - v*c];
        R2 = [1;v;0.5*v^2];
        R3 = [1;v + c;H + v*c];

        R = [R1,R2,R3];

        L = R\eye(3);
%         R = eye(3); L = eye(3);
        
        polx = L*polx; polxx = L*polxx; polxxx = L*polxxx;
        polxL = L*polxL; polxxL = L*polxxL; polxxxL = L*polxxxL;
        polxR = L*polx; polxxR = L*polxxR; polxxxR = L*polxxxR;
        polxmod = zeros(NumEq,1);
        polxxmod = zeros(NumEq,1);
        polxxxmod = zeros(NumEq,1);
        for n = 1:NumEq
            % reconstruction the polynomial on Ii 
            beta1 = 4*polxL(n)^2 + 16/5*polxL(n)*polxxxL(n) + 208/3*polxxL(n)^2 + 12496/5*polxxxL(n)^2;
            beta2 = 4*polx(n)^2 + 16/5*polx(n)*polxxx(n) + 208/3*polxx(n)^2 + 12496/5*polxxx(n)^2;
            beta3 = 4*polxR(n)^2 + 16/5*polxR(n)*polxxxR(n) + 208/3*polxxR(n)^2 + 12496/5*polxxxR(n)^2;
            tau5 = abs(beta3 - beta1)^2;
            omega1 = gammaN*(1 + tau5/(beta1 + epsilon));
            omega2 = gammaC*(1 + tau5/(beta2 + epsilon));
            omega3 = gammaN*(1 + tau5/(beta3 + epsilon));
            S = omega1 + omega2 + omega3;
            omega1 = omega1/S; omega2 = omega2/S; omega3 = omega3/S;
            polxmod(n) = omega1*polxL(n) + omega2*polx(n) + omega3*polxR(n);
            polxxmod(n) = omega1*polxxL(n) + omega2*polxx(n) + omega3*polxxR(n);
            polxxxmod(n) = omega1*polxxxL(n) + omega2*polxxx(n) + omega3*polxxxR(n);
        end
        polxmod = R*polxmod; polxxmod = R*polxxmod; polxxxmod = R*polxxxmod;
        for n = 1:NumEq
            uhmod(i,2,n) = polxmod(n);
            uhmod(i,3,n) = polxxmod(n);
            uhmod(i,4,n) = polxxxmod(n);
        end
    end
    
end

uh = uhmod;

end

