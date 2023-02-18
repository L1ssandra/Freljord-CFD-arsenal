% drawall.m
figure(1);
mesh(xc,yc,Q1);
title('rho')
colormap(cool);

figure(2);
mesh(xc,yc,Q2);
title('rhou')
colormap(cool);

figure(3);
mesh(xc,yc,Q3);
title('rhov')
colormap(cool);

figure(4);
mesh(xc,yc,Q4);
title('rhow')
colormap(cool);

figure(5);
mesh(xc,yc,Q5);
title('E')
colormap(cool);

figure(6);
mesh(xc,yc,Q6);
title('B1')
colormap(cool);

figure(7);
mesh(xc,yc,Q7);
title('B2')
colormap(cool);

figure(8);
mesh(xc,yc,Q8);
title('B3')
colormap(cool);