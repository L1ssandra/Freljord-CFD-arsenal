% drawall.m
figure(1);
mesh(rc,zc,Q1);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q1)),max(max(Q1)) + 1e-6]);
title('Bx')
colormap(cool);

figure(2);
mesh(rc,zc,Q2);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q2)),max(max(Q2)) + 1e-6]);
title('By')
colormap(cool);