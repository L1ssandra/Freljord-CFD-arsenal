%draw_solution1.m
figure(1)
hold on
plot(Xc,uh(:,1,1),'k-','linewidth',1.1)
plot(Xc1,u1(:,1,1),'c.-');
axis([-2,2,min(uh(:,1)) - 0.2,max(uh(:,1)) + 0.2]);
legend('exact','num')

figure(2)
hold on
uhv1 = u1(:,1,2)./u1(:,1,1);
uhv = uh(:,1,2)./uh(:,1,1);
plot(Xc,uhv,'k-','linewidth',1.1)
plot(Xc1,uhv1,'c.-'); 
axis([-2,2,min(uhv1) - 0.2,max(uhv1) + 0.2]);
legend('exact','num')

figure(3)
hold on
uhp1 = (gamma - 1)*(u1(:,1,3) - 0.5*u1(:,1,2).^2./u1(:,1,1));
uhp = (gamma - 1)*(uh(:,1,3) - 0.5*uh(:,1,2).^2./uh(:,1,1));
plot(Xc,uhp,'k-','linewidth',1.1)
plot(Xc1,uhp1,'c.-'); 
axis([-2,2,min(uhp) - 0.2,max(uhp) + 0.2]);
legend('exact','num')