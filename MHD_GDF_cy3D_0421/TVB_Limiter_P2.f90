    subroutine TVB_Limiter_P2
    
    include 'com.txt'
    
    real a00,a10,a01,a20,a11,a02,a30,a12
    real b00,b10,b01,b20,b11,b02,b21,b03
    real rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij
    integer t1,t2,t3
    
    call set_bc
    
    uhmod = uh
    
    uh(:,:,:,:,6:7) = 0
    
    do i = 0,Nx1
        do j = 0,Ny1
            uGint3D = 0
            do kk = 0,Nphi
                do n = 6,7
                    do d = 1,dimPk
                        uGint3D(:,:,kk,n) = uGint3D(:,:,kk,n) + uhmod(i,j,kk,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            do kk = 0,Nphi
                do n = 6,7
                    do d = 1,dimPk1
                        do i1 = 1,NumGLP
                            do j1 = 1,NumGLP
                                uh(i,j,kk,d,n) = uh(i,j,kk,d,n) + 0.25*weight(i1)*weight(j1)*uGint3D(i1,j1,kk,n)/(ra + (i - 0.5)*hr + hr1*lambda(i1))*phiG(i1,j1,d)/mm(d)
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    
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
                DeltaUR1(:,1) = uh(i,j,kk,2,:) + (2d0/3d0)*uh(i,j,kk,4,:)
                DeltaUL1(:,1) = uh(i,j,kk,2,:) - (2d0/3d0)*uh(i,j,kk,4,:)
            
                DeltaUR = matmul(L,DeltaUR)
                DeltaUL = matmul(L,DeltaUL)
                DeltaU = matmul(L,DeltaUR1)
            
                direction = 1
            
                call minmod
            
                DeltaUR1mod = matmul(R,DeltaUmod)
            
                do d = 1,NumEq
                    if (abs(DeltaUR1mod(d,1) - DeltaUR1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
                
                DeltaU = matmul(L,DeltaUL1)
                
                call minmod
                
                DeltaUL1mod = matmul(R,DeltaUmod)
                
                do d = 1,NumEq
                    if (abs(DeltaUL1mod(d,1) - DeltaUL1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
            
                ! y-direction
                call eigenmatrix(R,L,rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij,0.0,1.0,0.0)
                DeltaUR(:,1) = uh(i,j + 1,kk,1,:) - uh(i,j,kk,1,:)
                DeltaUL(:,1) = uh(i,j,kk,1,:) - uh(i,j - 1,kk,1,:)
                DeltaUU1(:,1) = uh(i,j,kk,3,:) + (2d0/3d0)*uh(i,j,kk,6,:)
                DeltaUD1(:,1) = uh(i,j,kk,3,:) - (2d0/3d0)*uh(i,j,kk,6,:)
            
                DeltaUR = matmul(L,DeltaUR)
                DeltaUL = matmul(L,DeltaUL)
                DeltaU = matmul(L,DeltaUU1)
            
                direction = 2
            
                call minmod
            
                DeltaUU1mod = matmul(R,DeltaUmod)
            
                do d = 1,NumEq
                    if (abs(DeltaUU1mod(d,1) - DeltaUU1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
                
                DeltaU = matmul(L,DeltaUD1)
                
                call minmod
                
                DeltaUD1mod = matmul(R,DeltaUmod)
                
                do d = 1,NumEq
                    if (abs(DeltaUD1mod(d,1) - DeltaUD1(d,1)) > 1e-6) then
                        change(d) = 1
                    end if
                end do
            
                do d = 1,NumEq
                    if (change(d) == 1) then
                        uhmod(i,j,kk,7:dimPk,d) = 0
                        uhmod(i,j,kk,5,d) = 0
                        uhmod(i,j,kk,2,d) = 0.5*(DeltaUR1mod(d,1) + DeltaUL1mod(d,1))
                        uhmod(i,j,kk,3,d) = 0.5*(DeltaUU1mod(d,1) + DeltaUD1mod(d,1))
                        uhmod(i,j,kk,4,d) = (3d0/4d0)*(DeltaUR1mod(d,1) - DeltaUL1mod(d,1))
                        uhmod(i,j,kk,6,d) = (3d0/4d0)*(DeltaUU1mod(d,1) - DeltaUD1mod(d,1))
                    end if
                end do
                
            end do
        end do
    end do
    
    uh = uhmod
    
    ! Limiting the magnetic field
    call set_bc
    
    do i = 0,Nx1
        do j = 0,Ny1
            do kk = 0,Nphi
                a00 = uh(i,j,kk,1,6)
                a10 = uh(i,j,kk,2,6)
                a01 = uh(i,j,kk,3,6)
                a20 = uh(i,j,kk,4,6)
                a11 = uh(i,j,kk,5,6)
                a02 = uh(i,j,kk,6,6)
                a30 = uh(i,j,kk,7,6)
                a12 = uh(i,j,kk,9,6)
                
                aRM(i,j,kk,2) = a01 + a11
                aRM(i,j,kk,3) = a02 + a12
                
                aLM(i,j,kk,2) = a01 - a11
                aLM(i,j,kk,3) = a02 - a12
            end do
        end do
    end do
    
    do i = 0,Nx1
        do j = 0,Ny1
            do kk = 0,Nphi
                b00 = uh(i,j,kk,1,7)
                b10 = uh(i,j,kk,2,7)
                b01 = uh(i,j,kk,3,7)
                b20 = uh(i,j,kk,4,7)
                b11 = uh(i,j,kk,5,7)
                b02 = uh(i,j,kk,6,7)
                b21 = uh(i,j,kk,8,7)
                b03 = uh(i,j,kk,10,7)
            
                bUM(i,j,kk,2) = b10 + b11
                bUM(i,j,kk,3) = b20 + b21
            
                bDM(i,j,kk,2) = b10 - b11
                bDM(i,j,kk,3) = b20 - b21
            end do
        end do
    end do
    
    do i = 0,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                do d = 2,k + 1
                    call minmodB(Brmod(i,j,kk,d),Br(i,j,kk,d),aRM(i,j,kk,d)*(ra + i*hr),aLM(i + 1,j,kk,d)*(ra + i*hr),M,hz)
                end do
            end do
        end do
    end do
    
    do i = 1,Nr
        do j = 0,Nz
            do kk = 0,Nphi
                do d = 2,k + 1
                    call minmodB(Bzmod(i,j,kk,d),Bz(i,j,kk,d),bUM(i,j,kk,d)*Rc(i),bDM(i,j + 1,kk,d)*Rc(i),M,hr)
                end do
            end do
        end do
    end do
    
    Br(:,:,:,2:k + 1) = Brmod(:,:,:,2:k + 1)
    Bz(:,:,:,2:k + 1) = Bzmod(:,:,:,2:k + 1)
    
    end subroutine TVB_Limiter_P2