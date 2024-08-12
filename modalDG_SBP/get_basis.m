% get_basis.m
global phiG phixG mm phiGR phiGL Pmat Lmat DNmat Nfmat
phiG = zeros(NumGLP,dimPk);
phixG = zeros(NumGLP,dimPk);
phiGR = zeros(1,dimPk);
phiGL = zeros(1,dimPk);
mm = zeros(1,dimPk);

for i = 1:NumGLP
    phiG(i,1) = 1;
    phiG(i,2) = lambda(i);
    phiG(i,3) = lambda(i)^2 - 1/3;
end

phiGR(1) = 1;
phiGR(2) = 1;
phiGR(3) = 2/3;

phiGL(1) = 1;
phiGL(2) = -1;
phiGL(3) = 2/3;

mm(1) = 1;
mm(2) = 1/3;
mm(3) = 4/45;

W = diag(weight); Wf = diag([1,1]);
Vq = phiG; Vf = [phiGR;phiGL];
Mmat = Vq'*W*Vq; Minv = Mmat\eye(dimPk);
D = [0, 0, 0; 1, 0, 0; 0, 2, 0]';
Nfmat = diag([1,-1]);
Pmat = Minv*Vq'*W; 
Lmat = Minv*Vf'*Wf;
Dq = Vq*D*Pmat;
DNmat = [ Dq - 0.5*Vq*Lmat*Nfmat*Vf*Pmat, 0.5*Vq*Lmat*Nfmat;  -0.5*Nfmat*Vf*Pmat,  0.5*Nfmat ];
Qmat = diag([weight,1,1])*DNmat;

% % test
% % Pmat*Vq
% Btest = Qmat + Qmat';
% % Btest = W*Dq + Dq'*W - Pmat'*Vf'*Wf*Nfmat*Vf*Pmat;
% % Btest = Mmat*D + D'*Mmat - Vf'*Wf*Nfmat*Vf;
% for i = 1:7
%     for j = 1:7
%         if abs(Btest(i,j)) < 1e-10
%             Btest(i,j) = 0;
%         end
%     end
% end
% Btest
