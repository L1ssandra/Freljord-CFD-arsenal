% drawB.m
figure(6);
mesh(rc,zc,Q6);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q6)),max(max(Q6)) + 1e-6]);
title('Br')
colormap(cool);

figure(7);
mesh(rc,zc,Q7);
axis([Rc(1),Rc(end),Zc(1),Zc(end),min(min(Q7)),max(max(Q7)) + 1e-6]);
title('Bz')
colormap(cool);