    module com
    
    integer Nx,Ny,frameMAX,Nx1,Ny1
    real(8) pi
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,maxu
    integer bcL,bcR,bcU,bcD
    
    parameter(Nx = 100, Ny = 100, frameMAX = 60)
    
    parameter(pi = 4*atan(1.0d0))
    parameter(Nx1 = Nx + 1, Ny1 = Ny + 1)

    real(8) hx,hy,X(0:Nx),Y(0:Ny)
    real(8) uh(-1:Nx1,-1:Ny1),du(-1:Nx1,-1:Ny1)
    
    end module com
    
    !*******************************
    
    module init1
    
    contains
    
    function u0(x,y)
    
    real(8) u0,x,y
    u0 = sin(x)*sin(y)
    
    end function u0
    
    subroutine mesh
    
    use com
    
    xa = 0
    xb = 2*pi
    ya = 0
    yb = 2*pi

    bcL = 1
    bcR = 1
    bcU = 1
    bcD = 1

    tend = 2*pi
    
    end subroutine mesh
    
    end module init1
    
    !*******************************
    
    program main
    
    use com
    
    call init_data
    
    call Euler_Forward
    
    end program main

    !*******************************
    
    subroutine init_data
    
    use com
    
    use init1
    call mesh
    
    hx = (xb - xa)/Nx
    hy = (yb - ya)/Ny
    
    uh = 0
    
    open(unit = 1,file = 'X.txt')
    open(unit = 2,file = 'Y.txt')
    
    do i = 0,Nx
        X(i) = xa + i*hx
        write(1,*) X(i)
    end do
    
    do j = 0,Ny
        Y(j) = ya + j*hy
        write(2,*) Y(j)
    end do
    
    do i = 0,Nx
        do j = 0,Ny
            uh(i,j) = u0(X(i),Y(j))
        end do
    end do
    
    close(1)
    close(2)
    
    end subroutine init_data

    !*******************************

    subroutine set_bc
    
    use com
    
    if (bcR == 1) then
        uh(Nx1,:) = uh(1,:)
    end if
    
    if (bcL == 1) then
        uh(-1,:) = uh(Nx - 1,:)
    end if
    
    if (bcU == 1) then
        uh(:,Ny1) = uh(:,1)
    end if
    
    if (bcD == 1) then
        uh(:,-1) = uh(:,Ny - 1)
    end if
    
    end subroutine set_bc


    !*******************************   


    subroutine Lh
    
    use com
    
    du = 0
    
    call set_bc
    
    do i = 0,Nx
        do j = 0,Ny
            du(i,j) = -(uh(i,j) - uh(i - 1,j))/hx - (uh(i,j) - uh(i,j - 1))/hy
        end do
    end do
    
    end subroutine Lh


    !******************************* 


    subroutine Euler_Forward
    
    use com
    
    real(8) t1
    
    CFL = 0.2
    dt = CFL*hx
    t = 0
    t1 = tend/frameMAX
    i1 = 1
    
    open(unit = 1,file = 'uh.txt')
    open(unit = 2,file = 'T.txt')
    
    do j = 0,Ny
        do i = 0,Nx
            write(1,*) uh(i,j)
        end do
    end do
    write(2,*) t
    
    do while (t < tend)
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        call Lh
        
        uh = uh + dt*du
        
        call calculate_maxu
        
        print *,t,maxu
        
        if (t >= i1*t1) then
            do j = 0,Ny
                do i = 0,Nx
                    write(1,*) uh(i,j)
                end do
            end do
            write(2,*) t
            print *,"save the solution at t = ",t,i1*t1
            i1 = i1 + 1
        end if
        
    end do
    
    end subroutine Euler_Forward

    !******************************* 

    subroutine calculate_maxu
    
    use com
    
    maxu = 0
    
    do i = 0,Nx
        do j = 0,Ny
            
            if (abs(uh(i,j)) > maxu) then
                maxu = abs(uh(i,j))
            end if
            
        end do
    end do
    
    end subroutine calculate_maxu


