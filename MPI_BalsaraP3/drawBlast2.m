% drawBlast2.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP2 = QB1.^2 + QB2.^2;
QV2 = Qu.^2 + Qv.^2;
xc = rc; yc = zc;

figure(1)
contour(xc,yc,Q1,40);
colormap(cool);
colorbar
title('Density')

figure(2);
contour(xc,yc,QP,40);
colormap(cool);
colorbar
title('Pressure')

figure(3);
contour(xc,yc,QV2,40);
colormap(cool);
colorbar
title('ux^2 + uy^2')

figure(4);
contour(xc,yc,QBP2,40);
colormap(cool);
colorbar
title('Bx^2 + By^2')