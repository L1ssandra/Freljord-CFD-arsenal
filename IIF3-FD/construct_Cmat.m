function Cmat = construct_Cmat(uh)

global Nx hx
Cmat = zeros(Nx,Nx);

for i = 1:Nx
    i1 = i + 1;
    if i1 > Nx
        i1 = 1;
    end
    Cmat(i,i) = -2;
    Cmat(i,i1) = 1;
    Cmat(i1,i) = 1;
end

Cmat = Cmat/hx^2;

end