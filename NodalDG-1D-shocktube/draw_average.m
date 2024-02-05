% draw_average.m

uibar = zeros(Nx,NumEq);

% calculate the cell-average
for i = 1:Nx
    for n = 1:NumEq
        uibar(i,n) = 0.5*uh(i,:,n)*Mvec;
    end
end

uh1 = uibar(:,1);
uh2 = uibar(:,2);
uh3 = uibar(:,3);

uhv = uh2./uh1;
uhp = (gamma - 1)*(uh3 - 0.5*uh1.*uhv.^2);

uh1 = log(uh1)/log(10);
uhp = log(uhp)/log(10);

Xq = zeros(1,Nx*NumGLP);
for i = 1:Nx
    Xq((i - 1)*NumGLP + 1:i*NumGLP) = Xc(i) + hx1*lambda;
end
figure(1); hold on
plot(Xc,uh1,'b-','linewidth',1.1);
axis([Xc(1),Xc(end),min(uh1) - 0.1,max(uh1) + 0.1]);

figure(2); hold on
plot(Xc,uhv,'b-','linewidth',1.1);
axis([Xc(1),Xc(end),min(uhv) - 0.1,max(uhv) + 0.1]);

figure(3); hold on
plot(Xc,uhp,'b-','linewidth',1.1);
axis([Xc(1),Xc(end),min(uhp) - 0.1,max(uhp) + 0.1]);