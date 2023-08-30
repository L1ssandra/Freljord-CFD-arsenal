% drawcontour.m
figure(1);
p = pcolor(rc,zc,Q1);p.EdgeColor = 'none';colormap(jet);colorbar;
% contour(rc,zc,Q1,30);colormap(cool)

% figure(2);
% gamma = 1.4;
% p = pcolor(rc,zc,(gamma - 1)*(Q4 - 0.5*(Q2.^2 + Q3.^2)./Q1));p.EdgeColor = 'none';colormap(jet);colorbar;
% %contour(xc,yc,Q1,40);colormap(cool)
% title('pressure')