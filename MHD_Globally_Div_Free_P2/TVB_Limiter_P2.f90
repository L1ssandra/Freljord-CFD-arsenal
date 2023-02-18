    subroutine TVB_Limiter_P2
    
    include 'com.txt'
    
    real a00,a10,a01,a20,a11,a02,a30,a12
    real b00,b10,b01,b20,b11,b02,b21,b03
    real rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij
    integer t1,t2,t3
    
    !call SYSTEM_CLOCK(t1)
    
    !call KXRCF_Detector
    
    !call SYSTEM_CLOCK(t2)
    
    uhmod = uh
    !print *,sum(Is_Trouble_Cell),Nx*Ny
    
    do i = 1,Nx
        do j = 1,Ny
            
            !if (Is_Trouble_Cell(i,j) == 1) then
                
                change = 0
            
                rhoij = uh(i,j,1,1)
                u1ij = uh(i,j,1,2)/rhoij
                u2ij = uh(i,j,1,3)/rhoij
                u3ij = uh(i,j,1,4)/rhoij
                Eij = uh(i,j,1,5)
                B1ij = uh(i,j,1,6)
                B2ij = uh(i,j,1,7)
                B3ij = uh(i,j,1,8)
            
                ! x-direction
                call eigenmatrix(R,L,rhoij,u1ij,u2ij,u3ij,Eij,B1ij,B2ij,B3ij,1.0,0.0,0.0)
                DeltaUR(:,1) = uh(i + 1,j,1,:) - uh(i,j,1,:)
                DeltaUL(:,1) = uh(i,j,1,:) - uh(i - 1,j,1,:)
                DeltaUR1(:,1) = uh(i,j,2,:) + (2d0/3d0)*uh(i,j,4,:) !- (1d0/3d0)*uh(i,j,6,:)
                DeltaUL1(:,1) = uh(i,j,2,:) - (2d0/3d0)*uh(i,j,4,:) !+ (1d0/3d0)*uh(i,j,6,:)
            
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
                call eigenmatrix(R,L,rhoij,u1ij,u2ij,0.0,Eij,B1ij,B2ij,0.0,0.0,1.0,0.0)
                DeltaUR(:,1) = uh(i,j + 1,1,:) - uh(i,j,1,:)
                DeltaUL(:,1) = uh(i,j,1,:) - uh(i,j - 1,1,:)
                DeltaUU1(:,1) = uh(i,j,3,:) + (2d0/3d0)*uh(i,j,6,:)!- (1d0/3d0)*uh(i,j,4,:) + (2d0/3d0)*uh(i,j,6,:)
                DeltaUD1(:,1) = uh(i,j,3,:) - (2d0/3d0)*uh(i,j,6,:)!+ (1d0/3d0)*uh(i,j,4,:) - (2d0/3d0)*uh(i,j,6,:)
            
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
                        uhmod(i,j,7:dimPk,d) = 0
                        uhmod(i,j,5,d) = 0
                        uhmod(i,j,2,d) = 0.5*(DeltaUR1mod(d,1) + DeltaUL1mod(d,1))
                        uhmod(i,j,3,d) = 0.5*(DeltaUU1mod(d,1) + DeltaUD1mod(d,1))
                        uhmod(i,j,4,d) = (2d0/3d0)*(DeltaUR1mod(d,1) - DeltaUL1mod(d,1)) !+ 0.5*(DeltaUU1mod(d,1) - DeltaUD1mod(d,1))
                        uhmod(i,j,6,d) = (2d0/3d0)*(DeltaUU1mod(d,1) - DeltaUD1mod(d,1)) !0.5*(DeltaUR1mod(d,1) - DeltaUL1mod(d,1)) + DeltaUU1mod(d,1) - DeltaUD1mod(d,1)
                    end if
                end do
                
            !end if
            
        end do
    end do
    
    uh = uhmod
    
    call pp_Limiter
    
    ! Limiting the magnetic field
    call set_bc
    
    do i = 0,Nx1
        do j = 0,Ny1
            a00 = uh(i,j,1,6)
            a10 = uh(i,j,2,6)
            a01 = uh(i,j,3,6)
            a20 = uh(i,j,4,6)
            a11 = uh(i,j,5,6)
            a02 = uh(i,j,6,6)
            a30 = uh(i,j,7,6)
            a12 = uh(i,j,9,6)
            
            aRM(i,j,1) = a00 + a10 + (2d0/3d0)*a20 + (2d0/5d0)*a30
            aRM(i,j,2) = a01 + a11
            aRM(i,j,3) = a02 + a12
            
            aLM(i,j,1) = a00 - a10 + (2d0/3d0)*a20 - (2d0/5d0)*a30
            aLM(i,j,2) = a01 - a11
            aLM(i,j,3) = a02 - a12
        end do
    end do
    
    do i = 0,Nx1
        do j = 0,Ny1
            b00 = uh(i,j,1,7)
            b10 = uh(i,j,2,7)
            b01 = uh(i,j,3,7)
            b20 = uh(i,j,4,7)
            b11 = uh(i,j,5,7)
            b02 = uh(i,j,6,7)
            b21 = uh(i,j,8,7)
            b03 = uh(i,j,10,7)
            
            bUM(i,j,1) = b00 + b01 + (2d0/3d0)*b02 + (2d0/5d0)*b03
            bUM(i,j,2) = b10 + b11
            bUM(i,j,3) = b20 + b21
            
            bDM(i,j,1) = b00 - b01 + (2d0/3d0)*b02 - (2d0/5d0)*b03
            bDM(i,j,2) = b10 - b11
            bDM(i,j,3) = b20 - b21
        end do
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            do d = 2,k + 1
                call minmodB(Bxmod(i,j,d),Bx(i,j,d),aRM(i,j,d),aLM(i + 1,j,d),M,hy)
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 2,k + 1
                call minmodB(Bymod(i,j,d),By(i,j,d),bUM(i,j,d),bDM(i,j + 1,d),M,hx)
            end do
        end do
    end do
    
    Bx(:,:,2:k + 1) = Bxmod(:,:,2:k + 1)
    By(:,:,2:k + 1) = Bymod(:,:,2:k + 1)
    
    !call SYSTEM_CLOCK(t3)
    
    end subroutine TVB_Limiter_P2
    
    