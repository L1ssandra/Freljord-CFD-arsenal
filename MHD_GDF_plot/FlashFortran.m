%FlashFortran.m
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

Q1flash = load('Q1flash.txt');
Q2flash = load('Q2flash.txt');
Q3flash = load('Q3flash.txt');
Q4flash = load('Q4flash.txt');
Q5flash = load('Q5flash.txt');
Q6flash = load('Q6flash.txt');
Q7flash = load('Q7flash.txt');
Q8flash = load('Q8flash.txt');
T = load('T.txt');

frame = length(T);

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1hflash = reshape(Q1flash,Nx,Ny,dimPk,frame);
Q2hflash = reshape(Q2flash,Nx,Ny,dimPk,frame);
Q3hflash = reshape(Q3flash,Nx,Ny,dimPk,frame);
Q4hflash = reshape(Q4flash,Nx,Ny,dimPk,frame);
Q5hflash = reshape(Q5flash,Nx,Ny,dimPk,frame);
Q6hflash = reshape(Q6flash,Nx,Ny,dimPk,frame);
Q7hflash = reshape(Q7flash,Nx,Ny,dimPk,frame);
Q8hflash = reshape(Q8flash,Nx,Ny,dimPk,frame);

Q1flash = Q1hflash(:,:,1,:);
Q2flash = Q2hflash(:,:,1,:);
Q3flash = Q3hflash(:,:,1,:);
Q4flash = Q4hflash(:,:,1,:);
Q5flash = Q5hflash(:,:,1,:);
Q6flash = Q6hflash(:,:,1,:);
Q7flash = Q7hflash(:,:,1,:);
Q8flash = Q8hflash(:,:,1,:);