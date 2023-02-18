% quadloadRotor.m
Xc = load('Xc.txt');
Yc = load('Yc.txt');

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');
Q7 = load('Q7.txt');
Q8 = load('Q8.txt');

Nx = length(Xc);
Ny = length(Yc);

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1 = reshape(Q1,Ny,Nx)';
Q2 = reshape(Q2,Ny,Nx)';
Q3 = reshape(Q3,Ny,Nx)';
Q4 = reshape(Q4,Ny,Nx)';
Q5 = reshape(Q5,Ny,Nx)';
Q6 = reshape(Q6,Ny,Nx)';
Q7 = reshape(Q7,Ny,Nx)';
Q8 = reshape(Q8,Ny,Nx)';

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

Xc = xc(:,1);
Yc = yc(1,:)';

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
                Q1c(i,j) = Q1c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q1(2*(i - 1) + i1,2*(j - 1) + j1);
                Q2c(i,j) = Q2c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q2(2*(i - 1) + i1,2*(j - 1) + j1);
                Q3c(i,j) = Q3c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q3(2*(i - 1) + i1,2*(j - 1) + j1);
                Q4c(i,j) = Q4c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q4(2*(i - 1) + i1,2*(j - 1) + j1);
                Q5c(i,j) = Q5c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q5(2*(i - 1) + i1,2*(j - 1) + j1);
                Q6c(i,j) = Q6c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q6(2*(i - 1) + i1,2*(j - 1) + j1);
                QPc(i,j) = QPc(i,j) + hx1*hy1*weight(i1)*weight(j1)*QP(2*(i - 1) + i1,2*(j - 1) + j1);
                QBPc(i,j) = QBPc(i,j) + hx1*hy1*weight(i1)*weight(j1)*QBP(2*(i - 1) + i1,2*(j - 1) + j1);
                QMachc(i,j) = QMachc(i,j) + hx1*hy1*weight(i1)*weight(j1)*QMach(2*(i - 1) + i1,2*(j - 1) + j1);
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
contour(xc,yc,Q1,20);colormap(cool);
%mesh(xc,yc,Q1);colormap(cool);
title('Density')

figure(2);
contour(xc,yc,QP,20);colormap(cool);
title('Pressure')

figure(3);
contour(xc,yc,(2*QBP).^0.5,10);colormap(cool);
title('Magnetic pressure')

figure(4);
contour(xc,yc,QMach,30);colormap(cool);
title('Mach number')
