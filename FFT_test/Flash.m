% Flash.m

h = figure(1);
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

QF = uh;
FRAME = frame;
t0 = T(end)/(FRAME - 1);

for i = 1:FRAME
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    plot(X,QF(:,j),'b-')
    axis([X(1),X(end),min(min(QF)) - 0.1,max(max(QF)) + 0.1]);
    pause(0.0001);
end