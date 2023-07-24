% test.m
X = 0:pi/4:7*pi/4;
uh = sin(X);
[ucos,usin] = pointtocoeff(uh);