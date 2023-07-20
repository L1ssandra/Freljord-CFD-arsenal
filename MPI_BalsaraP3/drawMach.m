% drawMach.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

XcM = [];
YcM = [];
i1 = Nr; i2 = 1;
j1 = Nz; j2 = 1;
for i = 1:Nr
    if Rc(i) > 0.3 && Rc(i) < 0.7
        XcM = [XcM,Rc(i)];
        i1 = min([i,i1]);
        i2 = max([i,i2]);
    end
end

for j = 1:Nz
    if Zc(j) > 0.3 && Zc(j) < 0.7
        YcM = [YcM,Zc(j)];
        j1 = min([j,j1]);
        j2 = max([j,j2]);
    end
end

QMachM = QMach(i1:i2,j1:j2);
NxM = length(XcM);
NyM = length(YcM);

xcM = zeros(NxM,NyM);
ycM = zeros(NxM,NyM);

for j = 1:NyM
    xcM(:,j) = XcM;
end

for i = 1:NxM
    ycM(i,:) = YcM';
end

ZcM = YcM;
figure(1)
QMach1D = 0.5*(QMachM(round(11*NxM/40),:) + QMachM(round(11*NxM/40) + 1,:));
plot(ZcM,QMach1D,'k-','linewidth',1.3);
% plot(Yc,Bx(Nx/10,:),'r-','linewidth',1.3)
% title('Mach cut at r = 1.41')
axis([ZcM(1),ZcM(end),min(QMach1D) - 0.2,max(QMach1D) + 0.2])

figure(2)
QMach1D2 = 0.5*(QMachM(round(0.6125*NxM),:) + QMachM(round(0.6125*NxM) + 1,:));
plot(ZcM,QMach1D2,'k-','linewidth',1.3);
% plot(Yc,Bx(Nx/10,:),'r-','linewidth',1.3)
% title('Mach cut at r = 1.41')
axis([ZcM(1),ZcM(end),min(QMach1D2) - 0.2,max(QMach1D2) + 0.2])