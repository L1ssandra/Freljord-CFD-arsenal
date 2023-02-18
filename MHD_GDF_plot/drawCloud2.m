% drawCloud2.m
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
QBnorm = sqrt(QB1.^2 + QB2.^2 + QB3.^2);
lnQ1 = log(abs(Q1));


M = max(max(lnQ1));
m = min(min(lnQ1));
Q1gray = (lnQ1 - m)/(M - m);
figure(1)
p = pcolor(xc,yc,lnQ1);
p.EdgeColor = 'none';
title('ln(rho)')
colormap(nclCM(203,200))

MP = max(max(QP));
mP = min(min(QP));
QPgray = (QP - mP)/(MP - mP);
figure(2)
p = pcolor(xc,yc,QP);
p.EdgeColor = 'none';
title('Pressure')
colormap(nclCM(203,200))

MB1 = max(max(QBnorm));
mB1 = min(min(QBnorm));
QBnormgray = (QBnorm - mB1)/(MB1 - mB1);
figure(3)
p = pcolor(xc,yc,QBnorm);
p.EdgeColor = 'none';
title('||B||')
colormap(nclCM(203,200))