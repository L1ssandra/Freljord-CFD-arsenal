function uhmod = MR(uh,uhL,uhR)

global gammaC gammaN
epsilon = 1e-10;
uhmod = uh;
II = uh(2);
III = uh(3);
IV = uh(4);
% p01 = a0
% p11
K11 = 1/gammaC;
% omega
beta01 = min( [4*uhL(2)^2,4*uhR(2)^2] );
beta11 = 4*(K11*II)^2;
tau1 = (beta01 - beta11)^2;
omega01 = gammaN*(1 + tau1/(beta01 + epsilon));
omega11 = gammaC*(1 + tau1/(beta11 + epsilon));
S = omega01 + omega11;
%omega01 = omega01/S;
omega11 = omega11/S;
% p12
K12 = omega11*K11;
% p22
K22 = (1 - gammaN*K12)/gammaC;
M22 = 1/gammaC;
% omega
beta12 = 4*(K12*II)^2;
beta22 = 4*(K22*II)^2 + 208/3*(M22*III)^2;
tau2 = (beta12 - beta22)^2;
omega12 = gammaN*(1 + tau2/(beta12 + epsilon));
omega22 = gammaC*(1 + tau2/(beta22 + epsilon));
S = omega12 + omega22;
omega12 = omega12/S; omega22 = omega22/S;
% p23
K23 = K22*omega22 + K12*omega12;
M23 = M22*omega22;
% p33
K33 = (1 - gammaN*K23)/gammaC;
M33 = (1 - gammaN*M23)/gammaC;
N33 = 1/gammaC;
% omega
beta23 = 4*(K23*II)^2 + 208/3*(M23*III)^2;
beta33 = 4*(K33*II)^2 + 16/5*(K33*N33*II*IV) + 208/3*(M33*III)^2 + 12496/5*(N33*IV)^2;
tau3 = (beta23 - beta33)^2;
omega23 = gammaN*(1 + tau3/(beta23 + epsilon));
omega33 = gammaC*(1 + tau3/(beta33 + epsilon));
S = omega23 + omega33;
omega23 = omega23/S; omega33 = omega33/S;
uhmod(2) = (omega23*K23 + omega33*K33)*II;
uhmod(3) = (omega23*M23 + omega33*M33)*III;
uhmod(4) =               (omega33*N33)*IV;

end