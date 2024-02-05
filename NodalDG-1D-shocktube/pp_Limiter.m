function uh = pp_Limiter(uh)
global Nx NumGLP NumEq Mvec

epsilon = 1e-13;
uibar = zeros(Nx,NumEq);

% calculate the cell-average
for i = 1:Nx
    for n = 1:NumEq
        uibar(i,n) = 0.5*uh(i,:,n)*Mvec;
    end
end

% Limiting the density
for i = 1:Nx
    theta = 1;
    for i1 = 1:NumGLP
        thetaq = min(abs([ (uibar(i,1) - epsilon)/(uibar(i,1) - uh(i,i1,1)), 1 ]));
        if thetaq < theta
            theta = thetaq;
        end
    end
    uh(i,:,1) = uibar(i,1) + (uh(i,:,1) - uibar(i,1))*theta;
end

% Limiting the pressure
for i = 1:Nx
    theta = 1;
    for i1 = 1:NumGLP
        pq = pressure(uh(i,i1,1),uh(i,i1,2),uh(i,i1,3));
        if pq < epsilon
            pt = @(t) pressure((1 - t)*uibar(i,1) + t*uh(i,i1,1),(1 - t)*uibar(i,2) + t*uh(i,i1,2),(1 - t)*uibar(i,3) + t*uh(i,i1,3)) - epsilon;
            tq = bisect(pt,0,1,epsilon);
            if tq < theta
                theta = tq;
            end
        end
    end
    
    if theta < 1
        theta = 0.95*theta;
    end
    
    for n = 1:NumEq
        uh(i,:,n) = uibar(i,n) + (uh(i,:,n) - uibar(i,n))*theta;
    end
end

end