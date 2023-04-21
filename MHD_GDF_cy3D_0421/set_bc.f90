    subroutine set_bc
    
    include 'com.txt'
    
    if (bcR == 1) then
        uh(Nx1,:,:,:,:) = uh(1,:,:,:,:)
    else if (bcR == 2) then
        do j = 1,Ny
            do kk = 0,Nphi
                do n = 1,NumEq
                    call evenex(uh(Nx1,j,kk,:,n),uh(Nx,j,kk,:,n))
                end do
            end do
        end do
    else if (bcR == 4) then
        uh(Nx1,:,:,:,:) = uh(1,:,:,:,:)
        uh(Nx1,:,:,:,6:7) = uh(1,:,:,:,6:7)*rb/ra
    else if (bcR == 5) then
        uh(Nx1,:,:,:,:) = uh(1,:,:,:,:)
        uh(Nx1,:,:,1,8) = 1d0/rb
        uh(Nx1,:,:,2:dimPk,8) = 0
    end if
    
    if (bcL == 1) then
        uh(0,:,:,:,:) = uh(Nx,:,:,:,:)
    else if (bcL == 2) then
        do j = 1,Ny
            do kk = 0,Nphi
                do n = 1,NumEq
                    call evenex(uh(0,j,kk,:,n),uh(1,j,kk,:,n))
                end do
            end do
        end do
    else if (bcL == 3) then
        do j = 1,Ny
            do kk = 0,Nphi
                call evenex(uh(0,j,kk,:,1),uh(1,j,kk,:,1))
                call oddex(uh(0,j,kk,:,2),uh(1,j,kk,:,2))
                call evenex(uh(0,j,kk,:,3),uh(1,j,kk,:,3))
                call evenex(uh(0,j,kk,:,4),uh(1,j,kk,:,4))
                call evenex(uh(0,j,kk,:,5),uh(1,j,kk,:,5))
                call oddex(uh(0,j,kk,:,6),uh(1,Ny1 - j,kk,:,6))
                call evenex(uh(0,j,kk,:,7),uh(1,j,kk,:,7))
                call evenex(uh(0,j,kk,:,8),uh(1,j,kk,:,8))
            end do
        end do
    else if (bcL == 4) then
        uh(0,:,:,:,:) = uh(Nx,:,:,:,:)
        uh(0,:,:,:,6:7) = uh(Nx,:,:,:,6:7)*ra/rb
    else if (bcL == 5) then
        uh(0,:,:,:,:) = uh(Nx,:,:,:,:)
        uh(0,:,:,1,8) = 1d0/ra
        uh(0,:,:,2:dimPk,8) = 0
    end if
    
    if (bcU == 1) then
        uh(:,Ny1,:,:,:) = uh(:,1,:,:,:)
    else if (bcU == 2) then
        uh(:,Ny1,:,:,:) = uh(:,Ny,:,:,:)
    end if
    
    if (bcD == 1) then
        uh(:,0,:,:,:) = uh(:,Ny,:,:,:)
    else if (bcD == 2) then
        uh(:,0,:,:,:) = uh(:,1,:,:,:)
    else if (bcD == 4) then
        uh(:,0,:,:,:) = uh(:,1,:,:,:)
        do i = 1,Nx
            if (ra + i*hr <= 1.5) then
                uh(i,0,:,2:dimPk,:) = 0
                !uh(i,0,:,1,1) = U1(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,2) = U2(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,3) = U3(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,4) = U4(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,5) = U5(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,6) = U6(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,7) = U7(Rc(i),Zc(1) - hz,Phi(1))
                !uh(i,0,:,1,8) = U8(Rc(i),Zc(1) - hz,Phi(1))
            end if
        end do
    end if
    
    end subroutine set_bc
    
    
    subroutine evenex(a,b)
    
    real a(6),b(6)
    
    a(1) = b(1)
    a(2) = -b(2)
    a(3) = b(3)
    a(4) = b(4)
    a(5) = -b(5)
    a(6) = b(6)
    
    end subroutine evenex
    
    
    subroutine oddex(a,b)
    
    real a(6),b(6)
    
    a(1) = -b(1)
    a(2) = b(2)
    a(3) = -b(3)
    a(4) = -b(4)
    a(5) = b(5)
    a(6) = -b(6)
    
    end subroutine oddex