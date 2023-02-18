    subroutine real_solution
    
    include 'com.txt'
    
    include 'init1.txt'
    
    do i = 1,Nx
        do j = 1,Ny
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    ureal(i,j,i1,j1,1) = U1(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,2) = U2(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,3) = U3(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,4) = U4(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,5) = U5(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,6) = U6(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,7) = U7(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                    ureal(i,j,i1,j1,8) = U8(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend)
                end do
            end do
        end do
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                Ez1real(i,j,j1) = U6(xa + i*hx - tend,Yc(j) + hy1*lambda(j1) - tend)
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                Ez2real(i,j,i1) = U7(Xc(i) + hx1*lambda(i1) - tend,ya + j*hy - tend)
            end do
        end do
    end do
    
    end subroutine real_solution