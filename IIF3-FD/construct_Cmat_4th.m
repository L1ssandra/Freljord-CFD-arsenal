function Cmat = construct_Cmat_4th(uh)

global Nx hx
Cmat = zeros(Nx,Nx);

for i = 1:Nx
    i1 = i + 1;
    i2 = i + 2;
    if i1 == Nx + 1
        i1 = 1;
    end
    if i2 == Nx + 1
        i2 = 1;
    elseif i2 == Nx + 2
        i2 = 2;
    end
    Cmat(i,i) = -30;
    Cmat(i,i1) = 16;
    Cmat(i,i2) = -1;
    Cmat(i1,i) = 16;
    Cmat(i2,i) = -1;
end

Cmat = Cmat/(12*hx^2);

end