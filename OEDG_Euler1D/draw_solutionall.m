uhrhoplot = uh(:,1,1);
uhvplot = uh(:,1,2)./uh(:,1,1);
uhpplot = (gamma - 1)*(uh(:,1,3) - 0.5*uh(:,1,2).^2./uh(:,1,1));

% Leblanc
uhrhoplot = log(uhrhoplot); uhpplot = log(uhpplot);

figure(1)
plot(Xc,uhrhoplot,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhrhoplot) - 0.2,max(uhrhoplot) + 0.2]);

figure(2)
plot(Xc,uhvplot,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhvplot) - 0.2,max(uhvplot) + 0.2]);

figure(3)
plot(Xc,uhpplot,'b-','linewidth',1.3); 
axis([Xc(1),Xc(end),min(uhpplot) - 0.5,max(uhpplot) + 0.5]);