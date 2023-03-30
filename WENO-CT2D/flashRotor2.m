%flashRotor2.m

xa = 0;xb = 1;ya = 0;yb = 1; 
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



h = figure();
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');


[xc,yc] = meshgrid(Xc,Yc);
xlabel('X');
ylabel('Y');
grid on;
%s = [X(1),X(end),Y(1),Y(end),min(min(min(Q(:,:,:,4)))) - 0.1,max(max(max(Q(:,:,:,4)))) + 0.1];
s = [Xc(1),Xc(end),Yc(1),Yc(end),min(min(min(Q1))) - 0.01,max(max(max(Q1))) + 0.01];
axis(s);
TT = 100;
t0 = T(end)/TT;
for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %mesh(yc,xc,Q1(:,:,j));
    contourf(yc,xc,QMach(:,:,j),15);
    %axis(s);
    colormap(hot)
    title(T(j))
    pause(0.0001);
end
%mesh(yc,xc,Q1(:,:,j))
%axis(s)
%plot(T,DIV);
% contour(Xc,Yc,QF(:,:,j,1),25);

% figure(2)
% colormap(hot)
% contourf(yc,xc,Q1(:,:,j),15);
% 
% figure(3)
% colormap(hot)
% contourf(yc,xc,QBP(:,:,j),15);
% 
% figure(4)
% colormap(hot)
% contourf(yc,xc,QP(:,:,j),15);
% 
% figure(5)
% colormap(hot)
% contourf(yc,xc,QMach(:,:,j),15);