    subroutine real_solution
    
    include 'com.txt'
    
    include 'init1.txt'
    
    uh(:,:,:,:,2) = 0
    uh(:,:,:,:,4) = 0
    
    do i = 1,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                do d = 1,dimPk
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            !uh(i,j,kk,d,1) = uh(i,j,kk,d,1) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U1(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            uh(i,j,kk,d,2) = uh(i,j,kk,d,2) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U2(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            !uh(i,j,kk,d,3) = uh(i,j,kk,d,3) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U3(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            uh(i,j,kk,d,4) = uh(i,j,kk,d,4) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U4(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            !uh(i,j,kk,d,5) = uh(i,j,kk,d,5) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U5(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            !uh(i,j,kk,d,6) = uh(i,j,kk,d,6) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U6(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            !uh(i,j,kk,d,7) = uh(i,j,kk,d,7) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U7(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                            !uh(i,j,kk,d,8) = uh(i,j,kk,d,8) + 0.25*weight(i1)*weight(j1)*phiG(i1,j1,d)*U8(Rc(i) + hr1*lambda(i1) - tRK,Zc(j) + hz1*lambda(j1) - tRK,Phi(kk))/mm(d)
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    call set_bc
    
    end subroutine real_solution