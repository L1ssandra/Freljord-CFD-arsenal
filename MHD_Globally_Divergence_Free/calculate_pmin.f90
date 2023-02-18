    subroutine calculate_pmin
    
    include 'com.txt'
    
    pmin = 1000000
    
    uG = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    p1 = pressure(uG(i,j,i1,j1,1),uG(i,j,i1,j1,2),uG(i,j,i1,j1,3),uG(i,j,i1,j1,4),uG(i,j,i1,j1,5),uG(i,j,i1,j1,6),uG(i,j,i1,j1,7),uG(i,j,i1,j1,8),gamma)
                    if (p1 < pmin) then
                        pmin = p1
                    end if
                end do
            end do
        end do
    end do
    
    end subroutine calculate_pmin