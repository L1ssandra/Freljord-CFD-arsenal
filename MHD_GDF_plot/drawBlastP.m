%drawBlastP.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 1.4;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP2 = QB1.^2 + QB2.^2;
QV2 = Qu.^2 + Qv.^2;

figure(1)
p = pcolor(xc,yc,Q1);
colormap(nclCM(424,200));
p.EdgeColor = 'none';
colorbar
title('Density')

figure(2);
p = pcolor(xc,yc,QP);
colormap(nclCM(424,200));
p.EdgeColor = 'none';
colorbar
title('Pressure')

figure(3);
p = pcolor(xc,yc,QV2);
colormap(nclCM(424,200));
p.EdgeColor = 'none';
colorbar
title('ux^2 + uy^2')

figure(4);
p = pcolor(xc,yc,QBP2);
colormap(nclCM(424,200));
p.EdgeColor = 'none';
colorbar
title('Bx^2 + By^2')