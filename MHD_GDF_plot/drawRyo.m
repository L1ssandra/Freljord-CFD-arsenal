% drawRyo.m
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
plot(Xc,Q1(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Q1(:,3)) - 0.1,max(Q1(:,3)) + 0.1];
axis(s);
%mesh(xc,yc,Q1);colormap(cool);
title('Density')

figure(2);
plot(Xc,QP(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(QP(:,3)) - 0.1,max(QP(:,3)) + 0.1];
axis(s);
%mesh(xc,yc,Q1);colormap(cool);
title('Pressure')

figure(3);
plot(Xc,Q7(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Q7(:,3)) - 0.1,max(Q7(:,3)) + 0.1];
axis(s);
title('By')

figure(4);
plot(Xc,Qu(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Qu(:,3)) - 0.1,max(Qu(:,3)) + 0.1];
axis(s);
title('ux')

figure(5);
plot(Xc,Qv(:,3),'b.');colormap(cool);
s = [Xc(1),Xc(end),min(Qv(:,3)) - 0.1,max(Qv(:,3)) + 0.1];
axis(s);
title('uy')