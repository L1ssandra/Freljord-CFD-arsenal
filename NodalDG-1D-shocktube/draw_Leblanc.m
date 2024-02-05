uh1 = reshape(uh(:,:,1)',Nx*NumGLP,1);
uh2 = reshape(uh(:,:,2)',Nx*NumGLP,1);
uh3 = reshape(uh(:,:,3)',Nx*NumGLP,1);

uhv = uh2./uh1;
uhp = (gamma - 1)*(uh3 - 0.5*uh1.*uhv.^2);

uh1 = log(uh1)/log(10);
uhp = log(uhp)/log(10);

Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end
figure(1); hold on
plot(Xq,uh1,'b-','linewidth',1.1);
axis([Xq(1),Xq(end),min(uh1) - 0.1,max(uh1) + 0.1]);

figure(2); hold on
plot(Xq,uhv,'b-','linewidth',1.1);
axis([Xq(1),Xq(end),min(uhv) - 0.1,max(uhv) + 0.1]);

figure(3); hold on
plot(Xq,uhp,'b-','linewidth',1.1);
axis([Xq(1),Xq(end),min(uhp) - 0.1,max(uhp) + 0.1]);


% axis([4,4.5,min(uhv) - 0.1,max(uhv) + 0.1]);