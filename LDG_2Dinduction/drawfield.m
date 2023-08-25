% drawfield.m
figure(1);
[q] = streamslice(Rc,Zc,Q1,Q2,2,'noarrows');
axis([Rc(1),Rc(end),Zc(1),Zc(end)])
for i = 1:length(q)
    q(i).LineWidth = 1.2;
    q(i).Color = 'b';
end
% q.ShowArrowHead = 'off';
% q.Marker = '-';
% contour(rc,zc,Q1,30);colormap(cool)