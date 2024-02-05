uhv = reshape(uh',Nx*NumGLP,1);
Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end
figure(1); hold on
plot(X,u,'b-','linewidth',1.3);
plot(Xq,uhv,'r--','linewidth',1.1);
% axis([Xq(1),Xq(end),min(uhv) - 0.1,max(uhv) + 0.1]);
axis([4,4.5,min(uhv) - 0.1,max(uhv) + 0.1]);