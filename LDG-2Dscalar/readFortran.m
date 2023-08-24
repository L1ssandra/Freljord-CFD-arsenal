% readFortran.m

Rc = load('Xc.txt');
Zc = load('Yc.txt');
dimPk = 1;
% T = load('T.txt');
Nr = length(Rc);
Nz = length(Zc);

% frame = length(T);

Q1 = load('u0.txt');

rc = zeros(Nr,Nz);
zc = zeros(Nr,Nz);

for j = 1:Nz
    rc(:,j) = Rc;
end

for i = 1:Nr
    zc(i,:) = Zc';
end

Q1h = reshape(Q1,Nr,Nz,dimPk);

Q1 = Q1h(:,:,1);