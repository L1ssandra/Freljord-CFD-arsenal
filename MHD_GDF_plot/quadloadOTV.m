% quadloadOTV.m
load('HLL-OTV192-1.mat')
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
% QE = Q4;
% QB1 = Q5;
% QB2 = Q6;
Qw = Q4./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
QB3 = Q8;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

Xc = xc(:,1);
Yc = yc(1,:)';

% lambda(1) = -0.8611363115940525752239465;
% lambda(2) = -0.3399810435848562648026658;
% lambda(3) = 0.3399810435848562648026658;
% lambda(4) = 0.8611363115940525752239465;
% 
% weight(1) = 0.3478548451374538573730639;
% weight(2) = 0.6521451548625461426269361;
% weight(3) = 0.6521451548625461426269361;
% weight(4) = 0.3478548451374538573730639;

lambda(1) = -0.5773502691896257645091488;
lambda(2) = 0.5773502691896257645091488;
        
weight(1) = 1;
weight(2) = 1;

Nx = length(Xc)/2;
Ny = length(Yc)/2;

Xcc = zeros(1,Nx);
Ycc = zeros(1,Ny);
Q1c = zeros(Nx,Ny);
Q2c = zeros(Nx,Ny);
Q3c = zeros(Nx,Ny);
Q4c = zeros(Nx,Ny);
Q5c = zeros(Nx,Ny);
Q6c = zeros(Nx,Ny);
QPc = zeros(Nx,Ny);
QBPc = zeros(Nx,Ny);
QMachc = zeros(Nx,Ny);

for i = 1:Nx
    Xcc(i) = sum(Xc(2*(i - 1) + 1:2*i))/2;
end

for j = 1:Ny
    Ycc(j) = sum(Yc(2*(j - 1) + 1:2*j))/2;
end

hx1 = (Xcc(2) - Xcc(1))/2;
hy1 = (Ycc(2) - Ycc(1))/2;

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:2
            for j1 = 1:2
                Q1c(i,j) = Q1c(i,j) + weight(i1)*weight(j1)*Q1(2*(i - 1) + i1,2*(j - 1) + j1);
                Q2c(i,j) = Q2c(i,j) + weight(i1)*weight(j1)*Q2(2*(i - 1) + i1,2*(j - 1) + j1);
                Q3c(i,j) = Q3c(i,j) + weight(i1)*weight(j1)*Q3(2*(i - 1) + i1,2*(j - 1) + j1);
                Q4c(i,j) = Q4c(i,j) + weight(i1)*weight(j1)*Q4(2*(i - 1) + i1,2*(j - 1) + j1);
                Q5c(i,j) = Q5c(i,j) + weight(i1)*weight(j1)*Q5(2*(i - 1) + i1,2*(j - 1) + j1);
                Q6c(i,j) = Q6c(i,j) + weight(i1)*weight(j1)*Q6(2*(i - 1) + i1,2*(j - 1) + j1);
                QPc(i,j) = QPc(i,j) + weight(i1)*weight(j1)*QP(2*(i - 1) + i1,2*(j - 1) + j1);
                QBPc(i,j) = QBPc(i,j) + weight(i1)*weight(j1)*QBP(2*(i - 1) + i1,2*(j - 1) + j1);
                QMachc(i,j) = QMachc(i,j) + weight(i1)*weight(j1)*QMach(2*(i - 1) + i1,2*(j - 1) + j1);
            end
        end
    end
end

xcc = zeros(Nx,Ny);
ycc = zeros(Nx,Ny);

for j = 1:Ny
    xcc(:,j) = Xcc;
end

for i = 1:Nx
    ycc(i,:) = Ycc';
end

Q1 = Q1c;
Q2 = Q2c;
Q3 = Q3c;
Q4 = Q4c;
Q5 = Q5c;
Q6 = Q6c;
QP = QPc;
QBP = QBPc;
QMach = QMachc;

xc = xcc;
yc = ycc;

figure(1);
contour(xc,yc,Q1,15);colormap(cool);
title('Density')

figure(2);
contour(xc,yc,QP,15);colormap(cool);
title('Pressure')