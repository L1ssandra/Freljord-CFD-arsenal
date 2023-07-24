function [dusin,ducos] = LhFFT(usin,ucos)

global f sinkX coskX L1 Nx fuh uc
dusin = zeros(1,L1);
ducos = zeros(1,L1);
uh = zeros(1,Nx);

for i = 1:Nx
    uh(i) = uc;
    for d = 1:L1
        uh(i) = uh(i) + sinkX(d,i)*usin(d) + coskX(d,i)*ucos(d);
    end
end

for i = 1:Nx
    fuh(i) = f(uh(i));
end

% Fourier transform of f(uh)
[fucos,fusin] = pointtocoeff(fuh);

for d = 1:L1
    dusin(d) = -d*fucos(d);
    ducos(d) = d*fusin(d);
end
dusin = -dusin;
ducos = -ducos;

end