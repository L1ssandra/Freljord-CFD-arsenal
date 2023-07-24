function [dusin,ducos] = Lh(usin,ucos)

global f sinkX coskX L1 Nx fuh uc L
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
fusin = zeros(1,L1);
fucos = zeros(1,L1);

for i = 1:Nx
    for d = 1:L1
        fusin(d) = fusin(d) + sinkX(d,i)*fuh(i);
        fucos(d) = fucos(d) + coskX(d,i)*fuh(i);
    end
end
fusin = fusin/L;
fucos = fucos/L;

for d = 1:L1
    dusin(d) = -d*fucos(d);
    ducos(d) = d*fusin(d);
end
dusin = -dusin;
ducos = -ducos;

end