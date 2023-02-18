% drawOTV.m
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
QBP = 0.5*(QB1.^2 + QB2.^2);

% for i = 1:Nx
%     for j = 1:Ny
%         if Qrho(i,j) < 1.127 || Qrho(i,j) > 5.587
%             Qrho(i,j) = NaN;
%         end
%     end
% end

figure(1);
contour(xc,yc,Qrho,15);
%mesh(xc,yc,Q1);
title('Density')
colormap(cool);

% figure(2)
% contour(xc,yc,QP,15);
% %mesh(xc,yc,Q1);
% %title('rho')
% colormap(cool);

% draw Bx at x = pi
figure(2)
[Nx,Ny] = size(Q6);
% Bx = load('Bx.txt');
% Bx = reshape(Bx,Nx/5,Ny);
plot(Yc,(QB1(Nx/2,:) + QB1(Nx/2 + 1,:))/2,'r-','linewidth',1.3);
% plot(Yc,Bx(Nx/10,:),'r-','linewidth',1.3)
title('Bx cut at x = \pi')
axis([Yc(1),Yc(end),min(QB1(Nx/2,:)) - 0.1,max(QB1(Nx/2,:)) + 0.1])

