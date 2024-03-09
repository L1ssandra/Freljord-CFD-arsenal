function uh = OE(uh)

global Nx dimPk bcL bcR NumEq phiGR phiGL k
global phixGR phixxGR phixGL phixxGL dt hx hx1

s = 1;
uhb = zeros(Nx + 2,dimPk,NumEq);
uhb(2:end - 1,:,:) = uh;

uhR = zeros(Nx + 1,NumEq);
uhL = zeros(Nx + 1,NumEq);
uhxR = zeros(Nx + 1,NumEq);
uhxL = zeros(Nx + 1,NumEq);
uhxxR = zeros(Nx + 1,NumEq);
uhxxL = zeros(Nx + 1,NumEq);

ujumpRL = zeros(Nx + 1,NumEq);
uxjumpRL = zeros(Nx + 1,NumEq);
uxxjumpRL = zeros(Nx + 1,NumEq);

sigmapre = zeros(Nx,dimPk,NumEq);
sigma = zeros(Nx,dimPk);
avguh = zeros(NumEq);
deno = zeros(NumEq);
coeff = zeros(dimPk);
epsilon = 1e-15;


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

% calculate the jump
for i = 1:Nx + 1
    for n = 1:NumEq
        for d = 1:dimPk
            uhR(i,n) = uhR(i,n) + uhb(i,d,n)*phiGR(d);
            uhxR(i,n) = uhxR(i,n) + uhb(i,d,n)*phixGR(d);
            uhxxR(i,n) = uhxxR(i,n) + uhb(i,d,n)*phixxGR(d);
            
            uhL(i,n) = uhL(i,n) + uhb(i + 1,d,n)*phiGL(d);
            uhxL(i,n) = uhxL(i,n) + uhb(i + 1,d,n)*phixGL(d);
            uhxxL(i,n) = uhxxL(i,n) + uhb(i + 1,d,n)*phixxGL(d);
        end
    end
end

for i = 1:Nx + 1
    for n = 1:NumEq
        ujumpRL(i,n) = abs(uhR(i,n) - uhL(i,n));
        uxjumpRL(i,n) = abs(uhxR(i,n) - uhxL(i,n));
        uxxjumpRL(i,n) = abs(uhxxR(i,n) - uhxxL(i,n));
    end
end

% calculate max(uh - uhbar) and the coeff (2m + 1)*hx^m/(2k - 1)m!
for n = 1:NumEq
    avguh(n) = sum(uh(:,1,n))/Nx;
end
for n = 1:NumEq
    deno(n) = max(abs(uh(:,1,n) - avguh(n)));
end
for dd = 1:dimPk
    d = dd - 1;
    coeff(dd) = (2*d + 1)*hx^d/((2*k - 1)*factorial(d));
end

% calculate the damping term
for i = 1:Nx
    for n = 1:NumEq
        sigmapre(i,1,n) = coeff(1)*(ujumpRL(i,n) + ujumpRL(i + 1,n))/(2*deno(n) + epsilon);
        sigmapre(i,2,n) = coeff(2)*(uxjumpRL(i,n) + uxjumpRL(i + 1,n))/(2*deno(n) + epsilon);
        sigmapre(i,3,n) = coeff(3)*(uxxjumpRL(i,n) + uxxjumpRL(i + 1,n))/(2*deno(n) + epsilon);
    end
    for d = 1:dimPk
        sigma(i,d) = max(sigmapre(i,d,:));
    end
    for d = 2:dimPk
        sigma(i,d) = sigma(i,d) + sigma(i,d - 1);
    end
end

% reduce the high-order term
for i = 1:Nx
    [alpha1,alpha2] = wavespeed(uh(i,1,:));
    alpha = s*max([abs(alpha1),abs(alpha2)]);
    for d = 2:dimPk
        for n = 1:NumEq
            uh(i,d,n) = uh(i,d,n)*exp( -alpha*dt/hx*sigma(i,d) );
        end
    end
end

% figure(1);plot(sigma(:,1))%plot(ujumpRL(:,1))%
% figure(2);plot(sigma(:,2) - sigma(:,1))%plot(uxjumpRL(:,1)*hx1)%
% figure(3);plot(sigma(:,3) - sigma(:,2));error%plot(uxxjumpRL(:,1)*hx1^2);error%

end