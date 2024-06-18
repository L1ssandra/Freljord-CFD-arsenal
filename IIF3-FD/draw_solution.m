figure(1);hold on
plot(Xc,ureal,'k-','linewidth',1)
plot(Xc,uh,'bo','linewidth',1)
axis([Xc(1),Xc(end),min(uh) - 0.1,max(uh) + 0.1]);