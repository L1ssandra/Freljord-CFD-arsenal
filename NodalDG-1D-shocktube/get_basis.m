% get_basis.m
global phiG phixG mm phiGR phiGL Mmat D B S Minv Mvec
phiG = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
    
    phixG(i,1) = 0;
    phixG(i,2) = 1/hx1;
    phixG(i,3) = 2*lambda(i)/hx1;
end

Minv = zeros(NumGLP,NumGLP);
Mmat = zeros(NumGLP,NumGLP);
Mvec = zeros(NumGLP,1);
D = zeros(NumGLP,NumGLP);
B = zeros(NumGLP,NumGLP);
for i = 1:NumGLP
    Mmat(i,i) = weight(i);
    Minv(i,i) = 1/weight(i);
    Mvec(i) = weight(i);
end
V = zeros(NumGLP,NumGLP);
for i = 1:NumGLP
    for j = 1:NumGLP
        V(i,j) = lambda(i)^(NumGLP - j);
    end
end
for i = 1:NumGLP
    for j = 1:NumGLP
        y = zeros(NumGLP,1); y(j) = 1;
        Lj = V\y;
        dLj = polyder(Lj);
        X = zeros(NumGLP - 1,1); X(1) = lambda(i)^(NumGLP - 2);
        for d = 2:NumGLP - 1
            X(d) = X(d - 1)/lambda(i);
            if isnan(X(d)) == 1
                X(d) = 0;
                X(end) = 1;
            end
        end
        D(i,j) = dLj*X;
    end
end
B(1,1) = -1; B(NumGLP,NumGLP) = 1;
S = Mmat*D;

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;