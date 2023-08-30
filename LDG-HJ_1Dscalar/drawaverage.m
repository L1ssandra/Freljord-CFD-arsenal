% drawaverage.m
hold on
plot(Xc,uh(:,1),'rx-');
axis([Xc(1),Xc(end),min(uh(:,1)) - 0.1,max(uh(:,1)) + 0.1])