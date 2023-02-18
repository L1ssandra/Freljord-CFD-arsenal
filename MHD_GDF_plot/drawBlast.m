% drawBlast.m

Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 1.4;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP2 = QB1.^2 + QB2.^2;
QV2 = Qu.^2 + Qv.^2;
QPn = QP;
for i = 1:Nx
    for j = 1:Ny
        if QPn(i,j) > 0
            QPn(i,j) = 0;
        end
    end
end

figure(1);
contour(xc,yc,Q1,40);colormap(bone);
%mesh(xc,yc,Q1);colormap(cool);
title('Density')

figure(2);
contour(xc,yc,QP,40);colormap(bone);
title('Pressure')

figure(3);
contour(xc,yc,QV2,40);colormap(bone);
title('ux^2 + uy^2')

figure(4);
contour(xc,yc,QBP2,40);colormap(bone);
title('B1^2 + B2^2')

[Nx,Ny] = size(QP);
QPN = QP;
for i = 1:Nx
    for j = 1:Ny
        if QP(i,j) > 0
            QPN(i,j) = 0;
        end
    end
end

figure(5);
contour(xc,yc,QPN,40);colormap(cool);
title('negative pressure')
colorbar