% readFortran.m

XcG = load('Xc.txt');
YcG = load('Yc.txt');
lambda = load('lambda.txt');
weight = load('weight.txt');
NumGLP = length(lambda);
NxG = length(XcG);
NyG = length(YcG);
Nx = NxG/NumGLP;
Ny = NyG/NumGLP;
dimPk = 10;
Xc = zeros(1,Nx);
Yc = zeros(1,Ny);

for i = 1:Nx
    Xc(i) = (XcG((i - 1)*NumGLP + 1) + XcG(i*NumGLP))/2;
end

for j = 1:Ny
    Yc(j) = (YcG((j - 1)*NumGLP + 1) + YcG(j*NumGLP))/2;
end

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');
Q7 = load('Q7.txt');
Q8 = load('Q8.txt');

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1h = reshape(Q1,Nx,Ny,dimPk);
Q2h = reshape(Q2,Nx,Ny,dimPk);
Q3h = reshape(Q3,Nx,Ny,dimPk);
Q4h = reshape(Q4,Nx,Ny,dimPk);
Q5h = reshape(Q5,Nx,Ny,dimPk);
Q6h = reshape(Q6,Nx,Ny,dimPk);
Q7h = reshape(Q7,Nx,Ny,dimPk);
Q8h = reshape(Q8,Nx,Ny,dimPk);

%drawRotor
%drawSmoothVortex
%drawOTV
%drawall
%drawdiv
Q1 = Q1h(:,:,1);
Q2 = Q2h(:,:,1);
Q3 = Q3h(:,:,1);
Q4 = Q4h(:,:,1);
Q5 = Q5h(:,:,1);
Q6 = Q6h(:,:,1);
Q7 = Q7h(:,:,1);
Q8 = Q8h(:,:,1);
draw1