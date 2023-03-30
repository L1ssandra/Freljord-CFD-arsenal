%draw4.m

xa = 0;xb = 2*pi;ya = 0;yb = 2*pi; 
% 网格
h = (xb - xa)/N;
h1 = h/2;
X = zeros(N,1);
Y = zeros(N,1);
for i = 1:N + 1
    X(i) = xa + (i - 1)*h;
    Y(i) = ya + (i - 1)*h;
end

% 整格点
Xc = (X(1:end - 1) + X(2:end))/2;
Yc = (Y(1:end - 1) + Y(2:end))/2;


[xc,yc] = meshgrid(Xc,Yc);

figure(1)
[~,j] = min(abs(T - 0.5));
colormap(cool)
contour(yc,xc,Q1(:,:,j,1),15);

figure(2)
[~,j] = min(abs(T - 2));
colormap(cool)
contour(yc,xc,Q1(:,:,j,1),15);

figure(3)
[~,j3] = min(abs(T - 3));
colormap(cool)
contour(yc,xc,Q1(:,:,j3,1),15);

figure(4)
[~,j] = min(abs(T - 4));
colormap(cool)
contour(yc,xc,Q1(:,:,j,1),15);

% figure(6)
% s = [Xc(1),Xc(end),0,4];
% plot(Yc,(QP(:,60,j3) + QP(:,61,j3))/2,'k-');
% title('Thermal Pressure at y = 0.625π and t = 3')
% axis(s)
% 
% figure(7)
% colormap(cool)
% contour(yc,xc,QP(:,:,j3,1),15);
% hold on
% plot([0,2*pi],[0.625*pi,0.625*pi],'r-','linewidth',1)
% title('Thermal Pressure at t = 3')