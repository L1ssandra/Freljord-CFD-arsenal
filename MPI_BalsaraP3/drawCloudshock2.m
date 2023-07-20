%drawCloudshock2.m
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

figure(1);
%contour(rc,zc,Q1,15);colormap(cool);
p = pcolor(rc,zc,Q1);colormap(bone);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
% title('Density')

figure(2)
%contour(rc,zc,QP,15);colormap(cool)
p = pcolor(rc,zc,QP);colormap(bone);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);

figure(3)
%contour(rc,zc,QMach,30);colormap(cool)
p = pcolor(rc,zc,Q6);colormap(bone);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);

figure(4)
%contour(rc,zc,QBP,15);colormap(cool)
p = pcolor(rc,zc,Q7);colormap(bone);p.EdgeColor = 'none';colorbar;
%mesh(xc,yc,Q1);
%title('rho')
% colormap(cool);