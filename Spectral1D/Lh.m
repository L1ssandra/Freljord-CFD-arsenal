function [dusin,ducos] = Lh(usin,ucos)

global f sinkX coskX L Nx fuh Nx1 uc
dusin = zeros(1,L);
ducos = zeros(1,L);
uh = zeros(1,Nx1);

for i = 1:Nx1
    uh(i) = uc;
    for d = 1:L
        uh(i) = uh(i) + sinkX(d,i)*usin(d) + coskX(d,i)*ucos(d);
    end
end

for i = 1:Nx1
    fuh(i) = f(uh(i));
end

% Fourier transform of f(uh)
fusin = zeros(1,L);
fucos = zeros(1,L);

for i = 1:Nx
    for d = 1:L
        fusin(d) = fusin(d) + sinkX(d,i)*fuh(i);
        fucos(d) = fucos(d) + coskX(d,i)*fuh(i);
    end
end
fusin = fusin/L;
fucos = fucos/L;

for d = 1:L
    dusin(d) = -d*fucos(d);
    ducos(d) = d*fusin(d);
end
dusin = -dusin;
ducos = -ducos;

end