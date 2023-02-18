%average_to_GLpoint.m
hx1 = (Xc(2) - Xc(1))/2;
hy1 = (Yc(2) - Yc(1))/2;
NewGLP = 10;

NumGLP = NewGLP;

lambda = [0.1488743389816312108848260;
    0.4333953941292471907992659;
    0.6794095682990244062343274;
    0.8650633666889845107320967;
    0.9739065285171717200779640;
    -0.1488743389816312108848260;
    -0.4333953941292471907992659;
    -0.6794095682990244062343274;
    -0.8650633666889845107320967;
    -0.9739065285171717200779640];

get_basis

XGflash = zeros(1,Nx*NumGLP);
YGflash = zeros(1,Ny*NumGLP);

for i = 1:Nx
    XGflash(NumGLP*(i - 1) + 1:NumGLP*i) = Xc(i) + hx1*lambda;
end

for j = 1:Ny
    YGflash(NumGLP*(j - 1) + 1:NumGLP*j) = Yc(j) + hy1*lambda;
end

xcG = zeros(Nx*NumGLP,Ny*NumGLP);
ycG = zeros(Nx*NumGLP,Ny*NumGLP);

for j = 1:Ny*NumGLP
    xcG(:,j) = XGflash;
end

for i = 1:Nx*NumGLP
    ycG(i,:) = YGflash;
end

Q1G = zeros(Nx*NumGLP,Ny*NumGLP);
Q2G = zeros(Nx*NumGLP,Ny*NumGLP);
Q3G = zeros(Nx*NumGLP,Ny*NumGLP);
Q4G = zeros(Nx*NumGLP,Ny*NumGLP);
Q5G = zeros(Nx*NumGLP,Ny*NumGLP);
Q6G = zeros(Nx*NumGLP,Ny*NumGLP);
Q7G = zeros(Nx*NumGLP,Ny*NumGLP);
Q8G = zeros(Nx*NumGLP,Ny*NumGLP);

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:NumGLP
            for j1 = 1:NumGLP
                for d = 1:dimPk
                    Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q1G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q1h(i,j,d)*phiG(i1,j1,d);
                    Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q2G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q2h(i,j,d)*phiG(i1,j1,d);
                    Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q3G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q3h(i,j,d)*phiG(i1,j1,d);
                    Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q4G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q4h(i,j,d)*phiG(i1,j1,d);
                    Q5G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q5G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q5h(i,j,d)*phiG(i1,j1,d);
                    Q6G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q6G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q6h(i,j,d)*phiG(i1,j1,d);
                    Q7G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q7G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q7h(i,j,d)*phiG(i1,j1,d);
                    Q8G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) = Q8G((i - 1)*NumGLP + i1,(j - 1)*NumGLP + j1) + Q8h(i,j,d)*phiG(i1,j1,d);
                end
            end
        end
    end
end

Q1 = Q1G;
Q2 = Q2G;
Q3 = Q3G;
Q4 = Q4G;
Q5 = Q5G;
Q6 = Q6G;
Q7 = Q7G;
Q8 = Q8G;

Xc = XGflash;
Yc = YGflash;
xc = xcG;
yc = ycG;