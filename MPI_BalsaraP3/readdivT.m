% readdivT.m

divT = load('divT.txt');

tplot = 0:50/(length(divT) - 1):50;
plot(tplot,divT,'k-','linewidth',1)
axis([0,50,0.95*min(divT),1.05*max(divT)])
xlabel('t')
ylabel('||divB||')