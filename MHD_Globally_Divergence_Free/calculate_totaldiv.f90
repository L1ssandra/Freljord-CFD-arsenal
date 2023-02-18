    subroutine calculate_totaldiv
    
    include 'com.txt'
    
    uGdiv = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimPk
                uGdiv(i,j,:,:) = uGdiv(i,j,:,:) + uh(i,j,d,6)*phixG(:,:,d) + uh(i,j,d,7)*phiyG(:,:,d)
            end do
        end do
    end do
    
    ! Gauss integral
    totaldiv = 0
    do i = 1,Nx
        do j = 1,Ny
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    totaldiv = totaldiv + weight(i1)*weight(j1)*uGdiv(i,j,i1,j1)**2
                end do
            end do
        end do
    end do
    
    totaldiv = (totaldiv*hx1*hy1)**0.5
    
    end subroutine calculate_totaldiv