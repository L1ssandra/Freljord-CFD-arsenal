    subroutine calculate_L2error
    
    include 'com.txt'
    
    !include 'init1.txt'
    
    L2 = 0
    uG = 0
    
    ! The value of num solution on the GL points
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    
    uE = (uG - ureal)**2
    
    do n = 1,NumEq
        do i = 1,Nx
            do j = 1,Ny
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        L2(n) = L2(n) + weight(i1)*weight(j1)*uE(i,j,i1,j1,n)
                    end do
                end do
            end do
        end do
    end do
    
    L2 = (L2*hx1*hy1)**0.5d0
    
    print *,L2
    
    L2 = 0
    BxG = 0
    ByG = 0
    
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                BxG(i,j,:) = BxG(i,j,:) + Bx(i,j,d)*EzG(:,d)
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                ByG(i,j,:) = ByG(i,j,:) + By(i,j,d)*EzG(:,d)
            end do
        end do
    end do
    
    B1E = (BxG - Ez1real)**2
    B2E = (ByG - Ez2real)**2
    
    do i = 1,Nx
        do j = 1,Ny
            do i1 = 1,NumGLP
                L2(5) = L2(5) + weight(i1)*B1E(i,j,i1)
                L2(6) = L2(6) + weight(i1)*B2E(i,j,i1)
            end do
        end do
    end do
    
    L2(5) = (L2(5)*hx1*hy1)**0.5d0
    L2(6) = (L2(6)*hx1*hy1)**0.5d0
    
    print *," "
    print *,L2(5),L2(6)
    
    end subroutine calculate_L2error