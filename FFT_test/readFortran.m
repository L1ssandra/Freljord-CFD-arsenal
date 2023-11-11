% readFortran.m

X = load('X.txt');
T = load('T.txt');
Nx = length(X) - 1;
Nx1 = Nx + 1;
frame = length(T);

uh = load('uh.txt');

uh = reshape(uh,Nx1,frame);

X(end + 1) = 2*pi;
uh(end + 1,:) = uh(1,:);

Flash