    subroutine TVB_Limiter_P1
    
    include 'com.txt'
    
    real a00,a10,a01,a20,a11,a02,a30,a12
    real b00,b10,b01,b20,b11,b02,b21,b03
    real rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij
    integer t1,t2,t3
    
    !call SYSTEM_CLOCK(t1)
    
    !call KXRCF_Detector
    
    !call SYSTEM_CLOCK(t2)
    
    uhmod = uh
    do i = 1,Nx
        do j = 1,Ny
            uh(i,j,:,:,6:7) = uh(i,j,:,:,6:7)/Rc(i)
        end do
    end do
    !print *,sum(Is_Trouble_Cell),Nx*Ny
    
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                
                change = 0
            
                rhoij = uh(i,j,kk,1,1)
                u1ij = uh(i,j,kk,1,2)/rhoij
                u2ij = uh(i,j,kk,1,3)/rhoij
                u3ij = uh(i,j,kk,1,4)/rhoij
                Eij = uh(i,j,kk,1,5)
                B1ij = uh(i,j,kk,1,6)
                B2ij = uh(i,j,kk,1,7)
                B3ij = uh(i,j,kk,1,8)
            
                ! x-direction
                call eigenmatrix(R,L,rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij,1.0,0.0,0.0)
                DeltaUR(:,1) = uh(i + 1,j,kk,1,:) - uh(i,j,kk,1,:)
                DeltaUL(:,1) = uh(i,j,kk,1,:) - uh(i - 1,j,kk,1,:)
                DeltaU(:,1) = uh(i,j,kk,2,:)
                DeltaU1 = DeltaU
            
                DeltaUR = matmul(L,DeltaUR)
                DeltaUL = matmul(L,DeltaUL)
                DeltaU = matmul(L,DeltaU)
            
                direction = 1
            
                call minmod
            
                DeltaUxmod = matmul(R,DeltaUmod)
            
                do d = 1,NumEq
                    if (abs(DeltaUxmod(d,1) - DeltaU1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
            
                ! y-direction
                call eigenmatrix(R,L,rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij,0.0,1.0,0.0)
                DeltaUR(:,1) = uh(i,j + 1,kk,1,:) - uh(i,j,kk,1,:)
                DeltaUL(:,1) = uh(i,j,kk,1,:) - uh(i,j - 1,kk,1,:)
                DeltaU(:,1) = uh(i,j,kk,3,:)
                DeltaU1 = DeltaU
            
                DeltaUR = matmul(L,DeltaUR)
                DeltaUL = matmul(L,DeltaUL)
                DeltaU = matmul(L,DeltaU)
            
                direction = 2
            
                call minmod
            
                DeltaUymod = matmul(R,DeltaUmod)
            
                do d = 1,NumEq
                    if (abs(DeltaUymod(d,1) - DeltaU1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
            
                do d = 1,NumEq
                    if (change(d) == 1) then
                        uhmod(i,j,kk,4:dimPk,d) = 0
                        uhmod(i,j,kk,2,d) = DeltaUxmod(d,1)
                        uhmod(i,j,kk,3,d) = DeltaUymod(d,1)
                    end if
                end do
                
            end do
        end do
    end do
    
    uh = uhmod
    
    call set_bc
    
    do i = 0,Nx1
        do j = 0,Ny1
            do kk = 0,Nphi
                a00 = uh(i,j,kk,1,6)
                a10 = uh(i,j,kk,2,6)
                a01 = uh(i,j,kk,3,6)
                a20 = uh(i,j,kk,4,6)
                a11 = uh(i,j,kk,5,6)
            
                aRM(i,j,kk,1) = a00 + a10 + (2d0/3d0)*a20
                aRM(i,j,kk,2) = a01 + a11
            
                aLM(i,j,kk,1) = a00 - a10 + (2d0/3d0)*a20
                aLM(i,j,kk,2) = a01 - a11
            end do
        end do
    end do
    
    do i = 0,Nx1
        do j = 0,Ny1
            do kk = 0,Nphi
                b00 = uh(i,j,kk,1,7)
                b10 = uh(i,j,kk,2,7)
                b01 = uh(i,j,kk,3,7)
                b11 = uh(i,j,kk,5,7)
                b02 = uh(i,j,kk,6,7)
            
                bUM(i,j,kk,1) = b00 + b01 + (2d0/3d0)*b02
                bUM(i,j,kk,2) = b10 + b11
            
                bDM(i,j,kk,1) = b00 - b01 + (2d0/3d0)*b02
                bDM(i,j,kk,2) = b10 - b11
            end do
        end do
    end do
    
    do i = 0,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                do d = 2,k + 1
                    call minmodB(Brmod(i,j,kk,d),Br(i,j,kk,d),aRM(i,j,kk,d),aLM(i + 1,j,kk,d),M,hz)
                end do
            end do
        end do
    end do
    
    do i = 1,Nr
        do j = 0,Nz
            do kk = 0,Nphi
                do d = 2,k + 1
                    call minmodB(Bzmod(i,j,kk,d),Bz(i,j,kk,d),bUM(i,j,kk,d),bDM(i,j + 1,kk,d),M,hr)
                end do
            end do
        end do
    end do
    
    Br(:,:,:,2:k + 1) = Brmod(:,:,:,2:k + 1)
    Bz(:,:,:,2:k + 1) = Bzmod(:,:,:,2:k + 1)
    
    end subroutine TVB_Limiter_P1