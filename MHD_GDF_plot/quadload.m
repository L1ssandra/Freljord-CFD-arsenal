% quadload.m
load('HLLC-Rotor200-1.mat')
Xc = xc(:,1);
Yc = yc(1,:)';

lambda(1) = -0.8611363115940525752239465;
lambda(2) = -0.3399810435848562648026658;
lambda(3) = 0.3399810435848562648026658;
lambda(4) = 0.8611363115940525752239465;

weight(1) = 0.3478548451374538573730639;
weight(2) = 0.6521451548625461426269361;
weight(3) = 0.6521451548625461426269361;
weight(4) = 0.3478548451374538573730639;

Nx = length(Xc)/4;
Ny = length(Yc)/4;

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

for i = 1:Nx
    Xcc(i) = sum(Xc(4*(i - 1) + 1:4*i))/4;
end

for j = 1:Ny
    Ycc(j) = sum(Yc(4*(j - 1) + 1:4*j))/4;
end

hx1 = (Xcc(2) - Xcc(1))/2;
hy1 = (Ycc(2) - Ycc(1))/2;

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:4
            for j1 = 1:4
                Q1c(i,j) = Q1c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q1(4*(i - 1) + i1,4*(j - 1) + j1);
                Q2c(i,j) = Q2c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q2(4*(i - 1) + i1,4*(j - 1) + j1);
                Q3c(i,j) = Q3c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q3(4*(i - 1) + i1,4*(j - 1) + j1);
                Q4c(i,j) = Q4c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q4(4*(i - 1) + i1,4*(j - 1) + j1);
                Q5c(i,j) = Q5c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q5(4*(i - 1) + i1,4*(j - 1) + j1);
                Q6c(i,j) = Q6c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q6(4*(i - 1) + i1,4*(j - 1) + j1);
                %Q7c(i,j) = Q7c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q7(4*(i - 1) + i1,4*(j - 1) + j1);
                %Q8c(i,j) = Q8c(i,j) + hx1*hy1*weight(i1)*weight(j1)*Q8(4*(i - 1) + i1,4*(j - 1) + j1);
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
%Q7 = Q7c;
%Q8 = Q8c;

xc = xcc;
yc = ycc;

drawRotor