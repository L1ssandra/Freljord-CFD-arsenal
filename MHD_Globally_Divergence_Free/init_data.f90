    subroutine init_data
    
    include 'com.txt'
    
    ! 0: Euler Vortex
    ! 1: MHD Vortex
    ! 2: Smooth Alfvén Wave 还没调好
    ! 3: Orszag-Tang Vortex
    ! 4: Rotor
    ! 5: MHD Blast
    ! 6: Cloud Shock Interaction
    ! 6b: Cloud Shock Interaction 2
    ! 7: Magnetic Loop
    ! 8: Brio-Wu Shock Tube
    ! 9: Ryu-Jones Shock Tube
    
    include 'init3.txt'
    
    hx = (xb - xa)/Nx
    hy = (yb - ya)/Ny
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    
    do i = 1,Nx
        Xc(i) = xa + (i - 0.5)*hx
    end do
    
    do j = 1,Ny
        Yc(j) = ya + (j - 0.5)*hy
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    ureal(i,j,i1,j1,1) = U1(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,2) = U2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,3) = U3(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,4) = U4(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,5) = U5(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,6) = U6(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,7) = U7(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                    ureal(i,j,i1,j1,8) = U8(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1))
                end do
            end do
        end do
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                Ez1real(i,j,j1) = U6(xa + i*hx,Yc(j) + hy1*lambda(j1))
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                Ez2real(i,j,i1) = U7(Xc(i) + hx1*lambda(i1),ya + j*hy)
            end do
        end do
    end do
    
    end subroutine init_data