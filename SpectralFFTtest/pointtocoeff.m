function [ucos,usin] = pointtocoeff(uh)

L = length(uh)/2;

% calculate the value of cL
CL = fft(uh)'/L;

L1 = L - 1;
ucos = zeros(1,L1);
usin = zeros(1,L1);
% calculate the value of aL and bL

for d = 1:L1
    ucos(d) = real((CL(1 + d) + CL(end + 1 - d))/2);
    usin(d) = real((CL(1 + d) - CL(end + 1 - d))/2i);
end

end
