function uh = MRWENO_Limiter(uh)

global dimPk Nx bcL bcR phiGL phiGR hx Ck gammaC gammaN
uhb = [zeros(1,dimPk);uh;zeros(1,dimPk)];
uhR = zeros(1,Nx + 1); uhL = zeros(1,Nx + 1);
ujump = zeros(1,Nx + 1);

uhmod = uh;
epsilon = 1e-10;

% set bc
if bcL == 1
    uhb(1,:) = uh(end,:);
end

if bcR == 1
    uhb(end,:) = uh(1,:);
end

% KXRCF Detector: find the trouble cells
for i = 1:Nx + 1
    for d = 1:dimPk
        uhR(i) = uhR(i) + uhb(i,d)*phiGR(d);
        uhL(i) = uhL(i) + uhb(i + 1,d)*phiGL(d);
    end
end
% calculate the jump
for i = 1:Nx + 1
    ujump(i) = uhL(i) - uhR(i);
end
% calculate the indicator
is_trouble_cell = zeros(1,Nx);
Ind = zeros(1,Nx);
for i = 1:Nx
    Ind(i) = (abs(ujump(i + 1)) + abs(ujump(i)))/(hx^2*max([abs(uhL(i)),abs(uhR(i))]) + 1e-16);
    if Ind(i) > Ck
        is_trouble_cell(i) = 1;
        %fprintf('%d  %d\n',Ind(i),i)
    end
end

for i = 1:Nx
    if is_trouble_cell(i) == 1
        
%         II = uh(i,2);
%         III = uh(i,3);
%         IV = uh(i,4);
%         % p01 = a0
%         % p11
%         K11 = 1/gammaC;
%         % omega
%         beta01 = min( [4*uhb(i,2)^2,4*uhb(i + 2,2)^2] );
%         beta11 = 4*(K11*II)^2;
%         tau1 = (beta01 - beta11)^2;
%         omega01 = gammaN*(1 + tau1/(beta01 + epsilon));
%         omega11 = gammaC*(1 + tau1/(beta11 + epsilon));
%         S = omega01 + omega11;
%         %omega01 = omega01/S; 
%         omega11 = omega11/S;
%         % p12
%         K12 = omega11*K11;
%         % p22
%         K22 = (1 - gammaN*K12)/gammaC;
%         M22 = 1/gammaC;
%         % omega
%         beta12 = 4*(K12*II)^2;
%         beta22 = 4*(K22*II)^2 + 208/3*(M22*III)^2;
%         tau2 = (beta12 - beta22)^2;
%         omega12 = gammaN*(1 + tau2/(beta12 + epsilon));
%         omega22 = gammaC*(1 + tau2/(beta22 + epsilon));
%         S = omega12 + omega22;
%         omega12 = omega12/S; omega22 = omega22/S;
%         % p23
%         K23 = K22*omega22 + K12*omega12;
%         M23 = M22*omega22;
%         % p33
%         K33 = (1 - gammaN*K23)/gammaC;
%         M33 = (1 - gammaN*M23)/gammaC;
%         N33 = 1/gammaC;
%         % omega
%         beta23 = 4*(K23*II)^2 + 208/3*(M23*III)^2;
%         beta33 = 4*(K33*II)^2 + 16/5*(K33*N33*II*IV) + 208/3*(M33*III)^2 + 12496/5*(N33*IV)^2;
%         tau3 = (beta23 - beta33)^2;
%         omega23 = gammaN*(1 + tau3/(beta23 + epsilon));
%         omega33 = gammaC*(1 + tau3/(beta33 + epsilon));
%         S = omega23 + omega33;
%         omega23 = omega23/S; omega33 = omega33/S;
%         uhmod(i,2) = (omega23*K23 + omega33*K33)*II;
%         uhmod(i,3) = (omega23*M23 + omega33*M33)*III;
%         uhmod(i,4) =               (omega33*N33)*IV;
%     end
        uhmod(i,:) = MR(uh(i,:),uhb(i,:),uhb(i + 2,:));
    
end

uh = uhmod;

end

