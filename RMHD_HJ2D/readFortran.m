% readFortran.m

Rc = load('Xc.txt');
Zc = load('Yc.txt');
dimPk = 1;
% T = load('T.txt');
Nr = length(Rc);
Nz = length(Zc);

% frame = length(T);

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');
Q7 = load('Q7.txt');
Q8 = load('Q8.txt');

Q1s = load('Q1s.txt');
Q2s = load('Q2s.txt');
Q3s = load('Q3s.txt');
Q4s = load('Q4s.txt');
Q5s = load('Q5s.txt');
Q6s = load('Q6s.txt');
Q7s = load('Q7s.txt');
Q8s = load('Q8s.txt');

rc = zeros(Nr,Nz);
zc = zeros(Nr,Nz);

for j = 1:Nz
    rc(:,j) = Rc;
end

for i = 1:Nr
    zc(i,:) = Zc';
end

Q1h = reshape(Q1,Nr,Nz,dimPk);
Q2h = reshape(Q2,Nr,Nz,dimPk);
Q3h = reshape(Q3,Nr,Nz,dimPk);
Q4h = reshape(Q4,Nr,Nz,dimPk);
Q5h = reshape(Q5,Nr,Nz,dimPk);
Q6h = reshape(Q6,Nr,Nz,dimPk);
Q7h = reshape(Q7,Nr,Nz,dimPk);
Q8h = reshape(Q8,Nr,Nz,dimPk);

Q1sh = reshape(Q1s,Nr,Nz,dimPk);
Q2sh = reshape(Q2s,Nr,Nz,dimPk);
Q3sh = reshape(Q3s,Nr,Nz,dimPk);
Q4sh = reshape(Q4s,Nr,Nz,dimPk);
Q5sh = reshape(Q5s,Nr,Nz,dimPk);
Q6sh = reshape(Q6s,Nr,Nz,dimPk);
Q7sh = reshape(Q7s,Nr,Nz,dimPk);
Q8sh = reshape(Q8s,Nr,Nz,dimPk);

Q1 = Q1h(:,:,1);
Q2 = Q2h(:,:,1);
Q3 = Q3h(:,:,1);
Q4 = Q4h(:,:,1);
Q5 = Q5h(:,:,1);
Q6 = Q6h(:,:,1);
Q7 = Q7h(:,:,1);
Q8 = Q8h(:,:,1);

Q1s = Q1sh(:,:,1);
Q2s = Q2sh(:,:,1);
Q3s = Q3sh(:,:,1);
Q4s = Q4sh(:,:,1);
Q5s = Q5sh(:,:,1);
Q6s = Q6sh(:,:,1);
Q7s = Q7sh(:,:,1);
Q8s = Q8sh(:,:,1);