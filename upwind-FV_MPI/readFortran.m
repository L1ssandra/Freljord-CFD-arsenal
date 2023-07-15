% readFortran.m

X = load('Xc.txt');
Y = load('Yc.txt');
% T = load('T.txt');
Nx = length(X);
Ny = length(Y);
Nx1 = Nx + 1;
Ny1 = Ny + 1;
frame = length(T);

uh = load('uh.txt');

x = zeros(Nx,Ny);
y = zeros(Nx,Ny);

for j = 1:Nx
    x(:,j) = X;
end

for i = 1:Ny
    y(i,:) = Y';
end

% uh = reshape(uh,Nx,Ny,frame);
uh = reshape(uh,Nx,Ny);