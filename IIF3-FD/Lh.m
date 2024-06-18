function du = Lh(uh)

global Nx bcL bcR hx %Xc t dt

N = Nx; h = hx; Q = uh';

fQ = zeros(N,1);

for i = 1:N
    fQ(i) = f(Q(i));
end

% 周期边界
if bcL == 1 && bcR == 1
    fQL2 = [fQ(end - 1:end);fQ(1:end - 2)];
    fQL1 = [fQ(end);fQ(1:end - 1)];
    fQR1 = [fQ(2:end);fQ(1)];
    fQR2 = [fQ(3:end);fQ(1:2)];
    fQR3 = [fQ(4:end);fQ(1:3)];

    QL2 = [Q(end - 1:end);Q(1:end - 2)];
    QL1 = [Q(end);Q(1:end - 1)];
    QR1 = [Q(2:end);Q(1)];
    QR2 = [Q(3:end);Q(1:2)];
    QR3 = [Q(4:end);Q(1:3)];
end

% 延长边界
% fQL2 = [fQ(1:2,:);fQ(1:end - 2,:)];
% fQL1 = [fQ(1,:);fQ(1:end - 1,:)];
% fQR1 = [fQ(2:end,:);fQ(end,:)];
% fQR2 = [fQ(3:end,:);fQ(end - 1:end,:)];
% fQR3 = [fQ(4:end,:);fQ(end - 2:end,:)];
% 
% QL2 = [Q(1:2,:);Q(1:end - 2,:)];
% QL1 = [Q(1,:);Q(1:end - 1,:)];
% QR1 = [Q(2:end,:);Q(end,:)];
% QR2 = [Q(3:end,:);Q(end - 1:end,:)];
% QR3 = [Q(4:end,:);Q(end - 2:end,:)];

fhat = zeros(1,N + 1);

alpha = 1.5;

for i = 1:N
    
    fQP = zeros(1,5);
    fQN = zeros(1,5);
    QP = zeros(1,5);
    QN = zeros(1,5);
    % f+
    fQP(1) = fQL2(i);
    fQP(2) = fQL1(i);
    fQP(3) = fQ(i);
    fQP(4) = fQR1(i);
    fQP(5) = fQR2(i);
    % f-
    fQN(1) = fQR3(i);
    fQN(2) = fQR2(i);
    fQN(3) = fQR1(i);
    fQN(4) = fQ(i);
    fQN(5) = fQL1(i);
    
    QP(1) = QL2(i);
    QP(2) = QL1(i);
    QP(3) = Q(i);
    QP(4) = QR1(i);
    QP(5) = QR2(i);
    
    QN(1) = QR3(i);
    QN(2) = QR2(i);
    QN(3) = QR1(i);
    QN(4) = Q(i);
    QN(5) = QL1(i);
    
    for j = 1:5
        fQP(j) = fQP(j);
        fQN(j) = fQN(j);
        QP(j) = QP(j);
        QN(j) = QN(j);
        fQP(j) = 0.5*(fQP(j) + alpha.*QP(j));
        fQN(j) = 0.5*(fQN(j) - alpha.*QN(j));
    end
    
    % 按分量重构
    QRnew = WENO5(fQP(1),fQP(2),fQP(3),fQP(4),fQP(5));
    QLnew = WENO5(fQN(1),fQN(2),fQN(3),fQN(4),fQN(5));
    
    Qnew = (QRnew + QLnew);
    
    fhat(i + 1) = Qnew;
    
end

% 周期边界
if bcL == 1 && bcR == 1
    fhat(1) = fhat(end);
end

% 延长边界
% fhat(1,:) = fhat(2,:);

du = -(fhat(2:end) - fhat(1:end - 1))/h;

% for i = 1:Nx
%     du(i) = -cos(Xc(i) - (t - dt));
% end

end
