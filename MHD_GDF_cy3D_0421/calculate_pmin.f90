    subroutine calculate_pmin
    
    include 'com.txt'
    
    pmin = 1000000
    
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                uGint = 0
                do n = 1,NumEq
                    do d = 1,dimPk
                        uGint(:,:,n) = uGint(:,:,n) + uh(i,j,kk,d,n)*phiG(:,:,d)
                    end do
                end do
                uhGLL = 0
                do d = 1,dimPk
                    do n = 1,NumEq
                        uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,kk,d,n)*phiGLL(:,:,d,:)
                    end do
                end do
            
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        do d = 1,2
                            p1 = pressure(uhGLL(i1,j1,1,d),uhGLL(i1,j1,2,d),uhGLL(i1,j1,3,d),uhGLL(i1,j1,4,d),uhGLL(i1,j1,5,d),uhGLL(i1,j1,6,d),uhGLL(i1,j1,7,d),uhGLL(i1,j1,8,d),gamma)
                            if (p1 < pmin) then
                                pmin = p1
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    end subroutine calculate_pmin
    
    
    function pressure(rho,rhou,rhov,rhow,E,B1,B2,B3,gamma)
    
    real rho,rhou,rhov,rhow,E,B1,B2,B3,gamma
    real pressure
    
    pressure = (gamma - 1)*(E - 0.5*(rhou**2 + rhov**2 + rhow**2)/rho - 0.5*(B1**2 + B2**2 + B3**2))
    
    !if (pressure < 0) then
    !    print *,E,(rhou**2 + rhov**2 + rhow**2)/rho,(B1**2 + B2**2 + B3**2)
    !end if
    
    end