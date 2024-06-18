% flash.m

figure(1)
t0 = T(end)/frameMAX;
for i = 1:frameMAX + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
   plot(Xc,uhflash(j,:),'b-');
   axis([Xc(1),Xc(end),min(min(uhflash)) - 0.1,max(max(uhflash)) + 0.1]);
%    axis([Xc(1),Xc(end),min(min(uhflash(j,:))) - 0.1,max(max(uhflash(j,:))) + 0.1]);
   pause(0.001);
end