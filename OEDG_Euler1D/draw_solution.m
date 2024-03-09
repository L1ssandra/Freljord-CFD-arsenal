%draw_solution.m
figure(1)
plot(Xc,uh(:,1,1),'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uh(:,1,1)) - 0.2,max(uh(:,1,1)) + 0.2]);

figure(2)
uhv = uh(:,1,2)./uh(:,1,1);
plot(Xc,uhv,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhv) - 0.2,max(uhv) + 0.2]);

figure(3)
uhp = (gamma - 1)*(uh(:,1,3) - 0.5*uh(:,1,2).^2./uh(:,1,1));
plot(Xc,uhp,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhp) - 0.2,max(uhp) + 0.2]);