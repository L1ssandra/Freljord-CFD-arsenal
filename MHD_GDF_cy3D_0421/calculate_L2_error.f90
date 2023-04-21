    subroutine calculate_L2_error
    
    include 'com.txt'
    
    include 'init7.txt'
    
    L2 = 0
    tend = 0
    
    ! The value of num solution on the GL points
    do i = 1,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                uGint = 0
                do n = 1,NumEq
                    do d = 1,dimPk
                        uGint(:,:,n) = uGint(:,:,n) + uh(i,j,kk,d,n)*phiG(:,:,d)
                    end do
                end do
            
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        L2(1) = L2(1) + weight(i1)*weight(j1)*(uGint(i1,j1,1) - U1(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(2) = L2(2) + weight(i1)*weight(j1)*(uGint(i1,j1,2) - U2(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(3) = L2(3) + weight(i1)*weight(j1)*(uGint(i1,j1,3) - U3(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(4) = L2(4) + weight(i1)*weight(j1)*(uGint(i1,j1,4) - U4(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(5) = L2(5) + weight(i1)*weight(j1)*(uGint(i1,j1,5) - U5(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(6) = L2(6) + weight(i1)*weight(j1)*(uGint(i1,j1,6)/(Rc(i) + hr1*lambda(i1)) - B1(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(7) = L2(7) + weight(i1)*weight(j1)*(uGint(i1,j1,7)/(Rc(i) + hr1*lambda(i1)) - B2(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                        L2(8) = L2(8) + weight(i1)*weight(j1)*(uGint(i1,j1,8) - U8(Rc(i) + hr1*lambda(i1) - tend,Zc(j) + hz1*lambda(j1) - tend,Phi(kk)))**2
                    end do
                end do
            end do
        end do
    end do
    
    L2 = (0.25*L2/(Nx*Ny*Nphi1))**0.5d0
    
    print *,L2
    
    end subroutine calculate_L2_error