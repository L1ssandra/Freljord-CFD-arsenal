figure(1);
plot(Xc,Q1(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Q1(:,3)) - 0.1,max(Q1(:,3)) + 0.1];
axis(s);
%mesh(xc,yc,Q1);colormap(cool);
title('Density')

figure(2);
plot(Xc,Q7(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Q7(:,3)) - 0.1,max(Q7(:,3)) + 0.1];
axis(s);
title('By')