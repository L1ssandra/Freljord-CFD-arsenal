% compare_time.m
Nplot = 100;
L = 64;
tic;
Spectral1D;
t2 = toc;
tic;
Spectral1DFFT;
t3 = toc;
fprintf('Spectral     : %d\n',t2)
fprintf('Spectral FFT : %d\n',t3)
