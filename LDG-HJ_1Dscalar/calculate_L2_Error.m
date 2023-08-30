% calculate_L2_Error.m
uhG = zeros(Nx,NumGLP);

for i = 1:Nx
    for d = 1:dimPk
        uhG(i,:) = uhG(i,:) + uh(i,d)*phiG(:,d)';
    end
end

ureal = zeros(Nx,NumGLP);

% Linear
for i = 1:Nx
    for j = 1:NumGLP
        ureal(i,j) = urealf(Xc(i) + hx1*lambda(j),tend);
    end
end

% Burgers
% main1;
% ureal = reshape(u,NumGLP,Nx)';

% sin
% for i = 1:Nx
%     for j = 1:NumGLP
%         ureal(i,j) = sin(2*atan(exp(-tend)*tan((Xc(i) + hx1*lambda(j))/2)));
%     end
% end

uE = abs(uhG - ureal);
L2_Error = 0;
L8_Error = max(max(uE));

for i = 1:Nx
    for i1 = 1:NumGLP
        L2_Error = L2_Error + hx1*weight(i1)*uE(i,i1)^2;
    end
end
L2_Error = sqrt(L2_Error);