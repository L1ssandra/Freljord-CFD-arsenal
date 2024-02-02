uhv = reshape(uh',Nx*NumGLP,1);
Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end
plot(Xq,uhv,'b-','linewidth',1.1);
axis([Xq(1),Xq(end),min(uhv) - 0.1,max(uhv) + 0.1]);