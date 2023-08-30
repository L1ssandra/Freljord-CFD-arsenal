function du = Lh(uh)
global Nx NumGLP dimPk mm phiG phixG phiGR phiGL bcL bcR weight hx Xc lambda hx1 eta
uhb = [zeros(1,dimPk);uh;zeros(1,dimPk)]; qhb = zeros(Nx + 2,dimPk);
qhG = zeros(Nx,NumGLP);
p1h = zeros(Nx,dimPk); p2h = zeros(Nx,dimPk);
uhG = zeros(Nx,NumGLP); p1hG = zeros(Nx,NumGLP); p2hG = zeros(Nx,NumGLP);
uhR = zeros(Nx + 1,1); qhR = zeros(Nx + 1,1);
uhL = zeros(Nx + 1,1); qhL = zeros(Nx + 1,1);
fhat = zeros(Nx + 1,1);
u1hat = zeros(Nx + 1,1);
u2hat = zeros(Nx + 1,1);
qhat = zeros(Nx + 1,1);
Bhat = zeros(Nx + 1,1);
du = zeros(Nx,dimPk);
HhatG = zeros(Nx,NumGLP);

% set_bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end

for i = 1:Nx + 1
    for d = 1:dimPk
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);
    end
end

% LDG: calculate the degree of p1 and p2.
% Here, we take p1 = ux+ and p2 = ux-.
for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    
    u1hat(i) = uR;
    u2hat(i) = uL;
    Bhat(i) = 0.5*(uR + uL);
end

for i = 1:Nx
    for d = 2:dimPk
        for i1 = 1:NumGLP
            p1h(i,d) = p1h(i,d) - 0.5*weight(i1)*uhG(i,i1)*phixG(i1,d);
            p2h(i,d) = p2h(i,d) - 0.5*weight(i1)*uhG(i,i1)*phixG(i1,d);
            qhb(i + 1,d) = qhb(i + 1,d) - 0.5*weight(i1)*uhG(i,i1)*phixG(i1,d);
        end
    end
end

for i = 1:Nx
    for d = 1:dimPk
        p1h(i,d) = p1h(i,d) + (1/hx)*(phiGR(d)*u1hat(i + 1) - phiGL(d)*u1hat(i));
        p2h(i,d) = p2h(i,d) + (1/hx)*(phiGR(d)*u2hat(i + 1) - phiGL(d)*u2hat(i));
        qhb(i + 1,d) = qhb(i + 1,d) + (1/hx)*(phiGR(d)*Bhat(i + 1) - phiGL(d)*Bhat(i));
    end
end

for d = 1:dimPk
    p1h(:,d) = p1h(:,d)/mm(d);
    p2h(:,d) = p2h(:,d)/mm(d);
    qhb(:,d) = qhb(:,d)/mm(d);
end

% set_bc for qh
if bcL == 1
    qhb(1,:) = qhb(Nx + 1,:);
end

if bcR == 1
    qhb(Nx + 2,:) = qhb(2,:);
end

% calculate the flux at edge.
for i = 1:Nx + 1
    for d = 1:dimPk
        qhR(i) = qhR(i) + qhb(i,d)*phiGR(d);
        qhL(i) = qhL(i) + qhb(i + 1,d)*phiGL(d);
    end
end

for i = 1:Nx + 1
    uR = uhL(i);
    uL = uhR(i);
    alpha = 1;
    fhat(i) = 0.5*(f(uR) + f(uL) - alpha*(uR - uL));
    
    qR = qhL(i);
    qL = qhR(i);
    qhat(i) = 0.5*(qR + qL);
end

% DG
for i = 1:Nx
    for d = 1:dimPk
        p1hG(i,:) = p1hG(i,:) + p1h(i,d)*phiG(:,d)';
        p2hG(i,:) = p2hG(i,:) + p2h(i,d)*phiG(:,d)';
        qhG(i,:) = qhG(i,:) + qhb(i + 1,d)*phiG(:,d)';
    end
end

% calculate Hhat
for i = 1:Nx
    for i1 = 1:NumGLP
        alpha = max([abs(dH(p1hG(i,i1))),abs(dH(p2hG(i,i1)))]);
%         alpha = max([abs(dH1(uhG(i,i1),p1hG(i,i1))),abs(dH1(uhG(i,i1),p2hG(i,i1)))]);
%         alpha = max([abs(dH2(p1hG(i,i1),Xc(i) + hx1*lambda(i1))),abs(dH2(p2hG(i,i1),Xc(i) + hx1*lambda(i1)))]);
        HhatG(i,i1) = H( 0.5*(p1hG(i,i1) + p2hG(i,i1)) ) - 0.5*alpha*(p1hG(i,i1) - p2hG(i,i1));
%         HhatG(i,i1) = H1( uhG(i,i1),0.5*(p1hG(i,i1) + p2hG(i,i1)) ) - 0.5*alpha*(p1hG(i,i1) - p2hG(i,i1));
%         HhatG(i,i1) = H2( 0.5*(p1hG(i,i1) + p2hG(i,i1)),Xc(i) + hx1*lambda(i1) ) - 0.5*alpha*(p1hG(i,i1) - p2hG(i,i1));
    end
end

for i = 1:Nx
    for d = 1:dimPk
        for i1 = 1:NumGLP
            du(i,d) = du(i,d) - 0.5*weight(i1)*(HhatG(i,i1)*phiG(i1,d) + eta*qhG(i,i1)*phixG(i1,d));
        end
    end
end

for i = 1:Nx
    for d = 1:dimPk
        du(i,d) = du(i,d) + (1/hx)*eta*(phiGR(d)*qhat(i + 1) - phiGL(d)*qhat(i));
    end
end


for d = 1:dimPk
    du(:,d) = du(:,d)/mm(d);
end

end

