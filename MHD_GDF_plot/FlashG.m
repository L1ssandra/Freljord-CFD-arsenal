% FlashG.m

h = figure();
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

QF = Q1flashG;
FRAME = 100;
t0 = T(end)/FRAME;
for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    %contour(xc,yc,QF(:,:,1,j),15);
    p = pcolor(xcG,ycG,QF(:,:,j));
    p.EdgeColor = 'none';
    %colormap(jet)
    colormap(nclCM(203,100));
    pause(0.0001);
end