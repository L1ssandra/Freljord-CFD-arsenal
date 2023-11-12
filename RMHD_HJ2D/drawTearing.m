% drawTearing.m
Qrho = Q1;
Qu = Q2;
Qv = Q3;
QP = Q5;
QB1 = Q6;
QB2 = Q7;

Qrhos = Q1s;
Qus = Q2s;
Qvs = Q3s;
QPs = Q5s;
QB1s = Q6s;
QB2s = Q7s;

figure(1);
p = pcolor(rc,zc,Q1 - Q1s);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,Q1 - Q1s,10,'linewidth',0.5);colormap(jet);colorbar;
title('rho')

figure(2);
p = pcolor(rc,zc,QP - QPs);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,QP - QPs,10,'linewidth',0.5);colormap(jet);colorbar;
title('p')

figure(3);
p = pcolor(rc,zc,Qu - Qus);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,Qu - Qus,10,'linewidth',0.5);colormap(jet);colorbar;
title('ux')

figure(4);
p = pcolor(rc,zc,Qv - Qvs);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,Qv - Qvs,10,'linewidth',0.5);colormap(jet);colorbar;
title('uy')

figure(5);
p = pcolor(rc,zc,QB1 - QB1s);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,QB1 - QB1s,10,'linewidth',0.5);colormap(jet);colorbar;
title('Bx')

figure(6);
p = pcolor(rc,zc,QB2 - QB2s);colormap(jet);p.EdgeColor = 'none';colorbar;
% contourf(rc,zc,QB2 - QB2s,10,'linewidth',0.5);colormap(jet);colorbar;
title('By')