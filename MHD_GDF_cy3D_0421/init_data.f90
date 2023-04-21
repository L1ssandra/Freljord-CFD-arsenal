    subroutine init_data
    
    include 'com.txt'
    
    ! 0: smooth vortex
    ! 1: sin-3D
    ! 2: divergence test
    ! 3: Orszag-Tang vortex
    ! 4: Orszag-Tang vortex 2
    ! 5: Spherical explosion
    ! 6: Rotor
    ! 7: Force-Balance-axisymmetric
    ! 8: Force-Balance-non-axisymmetric
    
    include 'init3.txt'
    
    hr = (rb - ra)/Nr
    hz = (zb - za)/Nz
    hphi = 2*pi/Nphi1
    hr1 = 0.5*hr
    hz1 = 0.5*hz
    
    do i = 1,Nr
        Rc(i) = ra + (i - 0.5)*hr
    end do
    
    do j = 1,Nz
        Zc(j) = za + (j - 0.5)*hz
    end do
    
    do kk = 0,Nphi
        Phi(kk) = kk*hphi
    end do
    
    call get_basis
    
    uh = 0
    
    ! L2 Pro
    do i = 1,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh(i,j,kk,d,1) = uh(i,j,kk,d,1) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U1(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,2) = uh(i,j,kk,d,2) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U2(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,3) = uh(i,j,kk,d,3) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U3(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,4) = uh(i,j,kk,d,4) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U4(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,5) = uh(i,j,kk,d,5) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U5(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,6) = uh(i,j,kk,d,6) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U6(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,7) = uh(i,j,kk,d,7) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U7(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                            uh(i,j,kk,d,8) = uh(i,j,kk,d,8) + weight(i1)*weight(j1)*phiG(i1,j1,d)*U8(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    uh = 0.25*uh
    
    do d = 1,dimPk1
        uh(:,:,:,d,:) = uh(:,:,:,d,:)/mm(d)
    end do
    
    ! The magnetic field on the edge
    Br = 0
    Bz = 0
    
    do i = 0,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                do d = 1,k + 1
                    do j1 = 1,NumGLP
                        Br(i,j,kk,d) = Br(i,j,kk,d) + 0.5*weight(j1)*U6(ra + i*hr,Zc(j) + hz1*lambda(j1),Phi(kk))*EzG(j1,d)
                    end do
                end do
            end do
        end do
    end do
    
    do i = 1,Nr
        do j = 0,Nz
            do kk = 0,Nphi
                do d = 1,k + 1
                    do i1 = 1,NumGLP
                        Bz(i,j,kk,d) = Bz(i,j,kk,d) + 0.5*weight(i1)*U7(Rc(i) + hr1*lambda(i1),za + j*hz,Phi(kk))*EzG(i1,d)
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,k + 1
        Br(:,:,:,d) = Br(:,:,:,d)/mmE(d)
        Bz(:,:,:,d) = Bz(:,:,:,d)/mmE(d)
    end do
    
    if (shock == 1) then
        do i = 1,Nx
            do j = 1,Ny
                do kk = 0,Nphi
                    uh(i,j,kk,2:dimPk,:) = 0
                    uh(i,j,kk,1,1) = U1(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,2) = U2(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,3) = U3(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,4) = U4(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,5) = U5(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,6) = U6(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,7) = U7(Rc(i),Zc(j),Phi(kk))
                    uh(i,j,kk,1,8) = U8(Rc(i),Zc(j),Phi(kk))
                end do
            end do
        end do
        
        do i = 0,Nr
            do j = 1,Nz
                do kk = 0,Nphi
                    Br(i,j,kk,2:k + 1) = 0
                    Br(i,j,kk,1) = U6(ra + i*hr,Zc(j),Phi(kk))
                end do
            end do
        end do
        
        do i = 1,Nr
            do j = 1,Nz
                do kk = 0,Nphi
                    Bz(i,j,kk,2:k + 1) = 0
                    Bz(i,j,kk,1) = U7(Rc(i),za + j*hz,Phi(kk))
                end do
            end do
        end do
        
    end if
    
    end subroutine init_data