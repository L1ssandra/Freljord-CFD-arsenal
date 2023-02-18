% quadtoaverage.m
% Xc = load('Xc.txt');
% Yc = load('Yc.txt');
% 
% Q1 = load('Q1.txt');
% Q2 = load('Q2.txt');
% Q3 = load('Q3.txt');
% Q4 = load('Q4.txt');
% Q5 = load('Q5.txt');
% Q6 = load('Q6.txt');
% Q7 = load('Q7.txt');
% Q8 = load('Q8.txt');

% Nx = length(Xc);
% Ny = length(Yc);
% 
% Q1 = reshape(Q1,Ny,Nx)';
% Q2 = reshape(Q2,Ny,Nx)';
% Q3 = reshape(Q3,Ny,Nx)';
% Q4 = reshape(Q4,Ny,Nx)';
% Q5 = reshape(Q5,Ny,Nx)';
% Q6 = reshape(Q6,Ny,Nx)';
% Q7 = reshape(Q7,Ny,Nx)';
% Q8 = reshape(Q8,Ny,Nx)';

Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
Qw = Q4./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
QB3 = Q8;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2 + Qw.^2) - 0.5*(QB1.^2 + QB2.^2 + QB3.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2 + Qw.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2 + QB3.^2);

NumGLP = 5;

lambda(1) = -0.9061798459386639927976269;
lambda(2) = -0.5384693101056830910363144;
lambda(3) = 0;
lambda(4) = 0.5384693101056830910363144;
lambda(5) = 0.9061798459386639927976269;

weight(1) = 0.2369268850561890875142640;
weight(2) = 0.4786286704993664680412915;
weight(3) = 0.5688888888888888888888889;
weight(4) = 0.4786286704993664680412915;
weight(5) = 0.2369268850561890875142640;

Nx = length(Xc)/NumGLP;
Ny = length(Yc)/NumGLP;

Xcc = zeros(1,Nx);
Ycc = zeros(1,Ny);
Q1c = zeros(Nx,Ny);
Q2c = zeros(Nx,Ny);
Q3c = zeros(Nx,Ny);
Q4c = zeros(Nx,Ny);
Q5c = zeros(Nx,Ny);
Q6c = zeros(Nx,Ny);
Q7c = zeros(Nx,Ny);
Q8c = zeros(Nx,Ny);
QPc = zeros(Nx,Ny);
QBPc = zeros(Nx,Ny);
QMachc = zeros(Nx,Ny);

for i = 1:Nx
    Xcc(i) = sum(Xc(NumGLP*(i - 1) + 1:NumGLP*i))/NumGLP;
end

for j = 1:Ny
    Ycc(j) = sum(Yc(NumGLP*(j - 1) + 1:NumGLP*j))/NumGLP;
end

hx1 = (Xcc(2) - Xcc(1))/2;
hy1 = (Ycc(2) - Ycc(1))/2;

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for j1 = 1:NumGLP
                Q1c(i,j) = Q1c(i,j) + 0.25*weight(i1)*weight(j1)*Q1(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q2c(i,j) = Q2c(i,j) + 0.25*weight(i1)*weight(j1)*Q2(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q3c(i,j) = Q3c(i,j) + 0.25*weight(i1)*weight(j1)*Q3(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q4c(i,j) = Q4c(i,j) + 0.25*weight(i1)*weight(j1)*Q4(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q5c(i,j) = Q5c(i,j) + 0.25*weight(i1)*weight(j1)*Q5(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q6c(i,j) = Q6c(i,j) + 0.25*weight(i1)*weight(j1)*Q6(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q7c(i,j) = Q7c(i,j) + 0.25*weight(i1)*weight(j1)*Q7(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                Q8c(i,j) = Q8c(i,j) + 0.25*weight(i1)*weight(j1)*Q8(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                QPc(i,j) = QPc(i,j) + 0.25*weight(i1)*weight(j1)*QP(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                QBPc(i,j) = QBPc(i,j) + 0.25*weight(i1)*weight(j1)*QBP(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
                QMachc(i,j) = QMachc(i,j) + 0.25*weight(i1)*weight(j1)*QMach(NumGLP*(i - 1) + i1,NumGLP*(j - 1) + j1);
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
Q7 = Q7c;
Q8 = Q8c;
QP = QPc;
QBP = QBPc;
QMach = QMachc;

xc = xcc;
yc = ycc;
Xc = Xcc;
Yc = Ycc;
                