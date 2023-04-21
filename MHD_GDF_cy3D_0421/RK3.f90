    subroutine RK3
    
    include 'com.txt'
    
    t = 0
    count = 0
    
    call TVB_Limiter_P2
    
    call div_free_Balsara
    
    call set_bc
    
    call calculate_umax
    
    call calculate_totaldiv
    
    call calculate_pmin
        
    print *,t,"  ",umax,"  ",totaldiv
    
    open(unit = 1,file = 'Q1flash.txt')
    open(unit = 2,file = 'Q2flash.txt')
    open(unit = 3,file = 'Q3flash.txt')
    open(unit = 4,file = 'Q4flash.txt')
    open(unit = 5,file = 'Q5flash.txt')
    open(unit = 6,file = 'Q6flash.txt')
    open(unit = 7,file = 'Q7flash.txt')
    open(unit = 8,file = 'Q8flash.txt')
    open(unit = 9,file = 'T.txt')
    
    do d = 1,dimPk
        do j = 1,Ny
            do i = 1,Nx
                write(1,*) uh(i,j,1,d,1)
                write(2,*) uh(i,j,1,d,2)
                write(3,*) uh(i,j,1,d,3)
                write(4,*) uh(i,j,1,d,4)
                write(5,*) uh(i,j,1,d,5)
                write(6,*) uh(i,j,1,d,6)/Rc(i)
                write(7,*) uh(i,j,1,d,7)/Rc(i)
                write(8,*) uh(i,j,1,d,8)
            end do
        end do
    end do
    write(9,*) t
    
    tt = tend/frameMAX
    countf = 1
    
    do while (t < tend)
        
        call calculate_dt
        
        tRK = t
        
        if (t + dt >= tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        if (dt < 1e-7) then
            t = tend
        end if
        
        ! Stage 1
        call set_bc
    
        call Lh
        
        call LhEz
        
        uI = uh + dt*du
        
        BrI = Br + dt*dBr
        BzI = Bz + dt*dBz
        
        uh0 = uh
        Br0 = Br
        Bz0 = Bz
        
        uh = uI
        Br = BrI
        Bz = BzI
        
        call TVB_Limiter_P2
        
        call div_free_Balsara
        
        if (RKorder == 3) then
            
            !Stage 2
            tRK = tRK + dt
            
            call set_bc
            
            call Lh
            
            call LhEz
            
            uh = uI + dt*du
            Br = BrI + dt*dBr
            Bz = BzI + dt*dBz
        
            uII = (3d0/4d0)*uh0 + (1d0/4d0)*uh
            BrII = (3d0/4d0)*Br0 + (1d0/4d0)*Br
            BzII = (3d0/4d0)*Bz0 + (1d0/4d0)*Bz
        
            uh = uII
            Br = BrII
            Bz = BzII
        
            call TVB_Limiter_P2
        
            call div_free_Balsara
            
            !Stage 3
            tRK = tRK - 0.5*dt
            
            call set_bc
            
            call Lh
            
            call LhEz
            
            Br = BrII + dt*dBr
            Bz = BzII + dt*dBz
            uh = uII + dt*du
        
            uh = (1d0/3d0)*uh0 + (2d0/3d0)*uh
            Br = (1d0/3d0)*Br0 + (2d0/3d0)*Br
            Bz = (1d0/3d0)*Bz0 + (2d0/3d0)*Bz
        
            call TVB_Limiter_P2
        
            call div_free_Balsara
            
        end if
        
        count = count + 1
        
        call calculate_umax
        
        call calculate_totaldiv
        
        call calculate_pmin
        
        print *,t,"  ",umax,"  ",totaldiv
        
        ! save solution
        if ((t > tt*countf .or. t == tend) .and. flash == 1) then
            print *,"save the solution at time",t,tt*countf
            do d = 1,dimPk
                do j = 1,Ny
                    do i = 1,Nx
                        write(1,*) uh(i,j,1,d,1)
                        write(2,*) uh(i,j,1,d,2)
                        write(3,*) uh(i,j,1,d,3)
                        write(4,*) uh(i,j,1,d,4)
                        write(5,*) uh(i,j,1,d,5)
                        write(6,*) uh(i,j,1,d,6)/Rc(i)
                        write(7,*) uh(i,j,1,d,7)/Rc(i)
                        write(8,*) uh(i,j,1,d,8)
                    end do
                end do
            end do
            write(9,*) t
            countf = countf + 1
        end if
        
    end do
    
    end subroutine RK3