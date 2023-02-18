    subroutine calculate_totaldiv
    
    include 'com.txt'
    
    totaldiv = 0
    
    do i = 1,Nx
        do j = 1,Ny
            uGdiv = 0
            do d = 1,dimPk
                uGdiv = uGdiv + uh(i,j,d,6)*phixG(:,:,d) + uh(i,j,d,7)*phiyG(:,:,d)
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    totaldiv = totaldiv + weight(i1)*weight(j1)*uGdiv(i1,j1)**2
                end do
            end do
            
        end do
    end do
    
    totaldiv = (totaldiv*hx1*hy1)**0.5
    
    end subroutine calculate_totaldiv