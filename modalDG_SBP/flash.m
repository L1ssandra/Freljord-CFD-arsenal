% flash.m

Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end

figure(1)
t0 = T(end)/frameMAX;
for i = 1:frameMAX + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
   plot(Xq,uhflash(j,:),'b-');
   axis([Xq(1),Xq(end),min(min(uhflash(j,:))) - 0.1,max(max(uhflash(j,:))) + 0.1]);
   pause(0.001);
end