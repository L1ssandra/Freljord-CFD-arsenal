    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny
    parameter(N_process = 100)
    parameter(Nx0 = 1200, Ny0 = 1200)
    parameter(Nx_process = 10, Ny_process = 10)
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,pi
    parameter(pi = 4*atan(1d0))
    integer bcR,bcL,bcU,bcD
    integer Right,Left,Up,Down
    
    real(8) uh(0:Nx1,0:Ny1),du(0:Nx1,0:Ny1),hx,hy,Xc(Nx),Yc(Ny),Xc0(Nx0),Yc0(Ny0)
    real(8) uh0(Nx0,Ny0)
    
    integer myid,myid1,the_id,the_id2
    integer myidx,myidy,the_idx,the_idy
    integer numprocs, namelen, rc,ierr,status(MPI_STATUS_SIZE),myid0
    character * (MPI_MAX_PROCESSOR_NAME) processor_name
    
    end module com
    
    
    !**********************************************************************
    
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

    tend = 1
    
    end subroutine mesh
    
    end module init1
    
    !**********************************************************************
    
    program main
    
    use com
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    
    myid1 = myid + 1
    
    myidx = mod(myid1,Nx_process)
    if (myidx == 0) then
        myidx = Nx_process
    end if
    
    myidy = (myid1 - myidx)/Nx_process + 1
    
    ! find the neighbor
    if (myidx /= Nx_process) then
        Right = myid + 1
    else
        Right = -100
    end if
    
    if (myidx /= 1) then
        Left = myid - 1
    else
        Left = -100
    end if
    
    if (myidy /= Ny_process) then
        Up = myid + Nx_process
    else
        Up = -100
    end if
    
    if (myidy /= 1) then
        Down = myid - Nx_process
    else
        Down = -100
    end if
    
    print *,"process",myid1,"is alive,the index is",myidx,myidy
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (myid1 == 1) then
        print *,""
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    print *,"process",myid1,"'s Right is",Right + 1
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (myid1 == 1) then
        print *,""
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    print *,"process",myid1,"'s Left is",Left + 1
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (myid1 == 1) then
        print *,""
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    print *,"process",myid1,"'s Up is",Up + 1
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if (myid1 == 1) then
        print *,""
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    print *,"process",myid1,"'s Down is",Down + 1
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    call init_data
    
    call Euler_Forward
    
    call save_solution
    
    call MPI_FINALIZE(rc)
    
    end program main
    
    !**********************************************************************
    
    subroutine init_data
    
    use com
    
    use init1
    call mesh
    
    hx = (xb - xa)/Nx0
    hy = (yb - ya)/Ny0
    
    uh0 = 0
    uh = 0
    
    do i = 1,Nx0
        Xc0(i) = xa + (i - 0.5)*hx
    end do
    
    do j = 1,Ny0
        Yc0(j) = ya + (j - 0.5)*hy
    end do
    
    do i = 1,Nx0
        do j = 1,Ny0
            uh0(i,j) = u0(Xc0(i),Yc0(j))
        end do
    end do
    
    uh(1:Nx,1:Ny) = uh0((myidx - 1)*Nx + 1:myidx*Nx,(myidy - 1)*Ny + 1:myidy*Ny)
    
    end subroutine init_data
    
    !**********************************************************************
    
    subroutine save_solution
    
    use com
    
    uh0 = 0
    uh0((myidx - 1)*Nx + 1:myidx*Nx,(myidy - 1)*Ny + 1:myidy*Ny) = uh(1:Nx,1:Ny)
    
    do the_id = 2,N_process
        
        i = mod(the_id,Nx_process)
        if (i == 0) then
            i = Nx_process
        end if
        j = (the_id - i)/Nx_process + 1
        
        if (myid1 == the_id) then
            call MPI_SEND(uh(1:Nx,1:Ny),Nx*Ny,MPI_REAL8,0,the_id,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(uh0((i - 1)*Nx + 1:i*Nx,(j - 1)*Ny + 1:j*Ny),Nx*Ny,MPI_REAL8,the_id - 1,the_id,MPI_COMM_WORLD,status,ierr) 
        end if
        
    end do
            
    if (myid1 == 1) then
        
        open(unit = 1,file = 'uh.txt')
        open(unit = 2,file = 'Xc.txt')
        open(unit = 3,file = 'Yc.txt')
        
        do j = 1,Ny0
            do i = 1,Nx0
                write(1,*) uh0(i,j)
            end do
        end do
        
        do i = 1,Nx0
            write(2,*) Xc0(i)
        end do
        
        do j = 1,Ny0
            write(3,*) Yc0(j)
        end do
        
        close(1)
        close(2)
        close(3)
        
    end if
    
    end subroutine save_solution
    
    !**********************************************************************
    
    subroutine Euler_Forward
    
    use com
    
    integer t1,t2
    
    CFL = 0.2
    dt = CFL*hx
    t = 0
    
    call system_clock(t1)
    
    do while (t < tend)
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        call Lh
        
        uh = uh + dt*du
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    call system_clock(t2)
    
    if (myid1 == 1) then
        open(unit = 1,file = 'time.txt')
        write(1,*) (t2 - t1)/10000.
        print *,(t2 - t1)/10000.
    end if
    
    end subroutine Euler_Forward
    
    !**********************************************************************
    
    subroutine Lh
    
    use com
    
    call set_bc1
    
    do i = 1,Nx
        do j = 1,Ny
            du(i,j) = -(uh(i,j) - uh(i - 1,j))/hx - (uh(i,j) - uh(i,j - 1))/hy
        end do
    end do
    
    end subroutine Lh
    
    !**********************************************************************
    
    subroutine set_bc
    
    use com
    
    do i = 1,Nx_process
        do j = 1,Ny_process
            
            the_id = i + Nx_process*(j - 1)
            
            ! The Right condition
            if (i == Nx_process) then
                the_idx = 1
                the_idy = j
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (bcR == 1) then
                    if (myid1 == the_id2) then
                        call MPI_SEND(uh(1,1:Ny),Ny,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                    end if
                end if
            
            else
                
                the_idx = i + 1
                the_idy = j
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (myid1 == the_id2) then
                    call MPI_SEND(uh(1,1:Ny),Ny,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                end if
                
            end if
            
            if (myid1 == the_id) then
                call MPI_RECV(uh(Nx1,1:Ny),Ny,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
            end if
            
            ! The Left condition
            if (i == 1) then
                the_idx = Nx_process
                the_idy = j
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (bcL == 1) then
                    if (myid1 == the_id2) then
                        call MPI_SEND(uh(Nx,1:Ny),Ny,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                    end if
                end if
            else
                the_idx = i - 1
                the_idy = j
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (myid1 == the_id2) then
                    call MPI_SEND(uh(Nx,1:Ny),Ny,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                end if
            end if
            
            if (myid1 == the_id) then
                call MPI_RECV(uh(0,1:Ny),Ny,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
            end if
            
            ! The Up condition
            if (j == Ny_process) then
                the_idx = i
                the_idy = 1
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (bcU == 1) then
                    if (myid1 == the_id2) then
                        call MPI_SEND(uh(1:Nx,1),Nx,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                    end if
                end if
            else
                the_idx = i
                the_idy = j + 1
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (myid1 == the_id2) then
                    call MPI_SEND(uh(1:Nx,1),Nx,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                end if
                
            end if
            
            if (myid1 == the_id) then
                call MPI_RECV(uh(1:Nx,Ny1),Nx,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
            end if
            
            ! The Down condition
            if (j == 1) then
                the_idx = i
                the_idy = Ny_process
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (bcD == 1) then
                    if (myid1 == the_id2) then
                        call MPI_SEND(uh(1:Nx,Ny),Nx,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                    end if
                end if
            else
                the_idx = i
                the_idy = j - 1
                the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                if (myid1 == the_id2) then
                    call MPI_SEND(uh(1:Nx,Ny),Nx,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                end if
                
            end if
            
            if (myid1 == the_id) then
                call MPI_RECV(uh(1:Nx,0),Nx,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
            end if
            
        end do
    end do
    
    end subroutine set_bc
    
    !**********************************************************************
    
    subroutine calculate_umax
    
    use com
    
    umax = 0
    umax1 = 0
    
    do i = 1,Nx
        do j = 1,Ny
            if (abs(uh(i,j)) > umax) then
                umax = abs(uh(i,j))
            end if
        end do
    end do
    
    do the_id = 2,N_process
        
        if (myid1 == the_id) then
            call MPI_SEND(umax,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(umax1,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            if (umax1 > umax) then
                umax = umax1
            end if
        end if
        
    end do
    
    do the_id = 2,N_process
        
        if (myid1 == 1) then
            call MPI_SEND(umax,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == the_id) then
            call MPI_RECV(umax,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
        end if
        
    end do
    
    end subroutine calculate_umax
    
    !**********************************************************************
    
    subroutine set_bc1
    
    use com
    
    use init1
    
    real(8) uhsendR(Ny),uhsendL(Ny),uhrecvR(Ny),uhrecvL(Ny)
    real(8) uhsendU(Nx),uhsendD(Nx),uhrecvU(Nx),uhrecvD(Nx)
    
    ! set bc
    if (myidx == 1) then ! Left
        do j = 1,Ny
            uh(0,j) = u0(Xc0(1) - hx - t,Yc0((myidy - 1)*Ny_process + j) - t)
        end do
    end if
    
    if (myidx == Nx_process) then ! Right
        do j = 1,Ny
            uh(Nx1,j) = u0(Xc0(Nx0) + hx - t,Yc0((myidy - 1)*Ny_process + j) - t)
        end do
    end if
    
    if (myidy == 1) then ! Down
        do i = 1,Nx
            uh(i,0) = u0(Xc0((myidx - 1)*Nx_process + i) - t,Yc0(0) - hy - t)
        end do
    end if
    
    if (mydiy == Ny_process) then ! Up
        do j = 1,Ny
            uh(i,Ny1) = u0(Xc0((myidx - 1)*Nx_process + i) - t,Yc0(Ny1) + hy - t)
        end do
    end if
    
    ! exchange bc
    ! odd  process: SEND first ---> RECV next
    ! even process: RECV first ---> SEND next
    
    ! x-direction
    if (mod(myidx,2) == 1) then ! odd
        
        if (Right >= 0) then
            uhsendR = uh(Nx,1:Ny)
            call MPI_SEND(uhsendR,Ny,MPI_REAL8,Right,myid,MPI_COMM_WORLD,ierr)
            call MPI_RECV(uhrecvR,Ny,MPI_REAL8,Right,Right + Nx0*Ny0,MPI_COMM_WORLD,status,ierr)
            uh(Nx1,1:Ny) = uhrecvR
        end if
    
        if (Left >= 0) then
            uhsendL = uh(1,1:Ny)
            call MPI_SEND(uhsendL,Ny,MPI_REAL8,Left,myid + Nx0*Ny0,MPI_COMM_WORLD,ierr)
            call MPI_RECV(uhrecvL,Ny,MPI_REAL8,Left,Left,MPI_COMM_WORLD,status,ierr)
            uh(0,1:Ny) = uhrecvL
        end if
        
    else ! even
        
        if (Right >= 0) then
            call MPI_RECV(uhrecvR,Ny,MPI_REAL8,Right,Right + Nx0*Ny0,MPI_COMM_WORLD,status,ierr)
            uh(Nx1,1:Ny) = uhrecvR
            uhsendR = uh(Nx,1:Ny)
            call MPI_SEND(uhsendR,Ny,MPI_REAL8,Right,myid,MPI_COMM_WORLD,ierr)
        end if
        
        if (Left >= 0) then
            call MPI_RECV(uhrecvL,Ny,MPI_REAL8,Left,Left,MPI_COMM_WORLD,status,ierr)
            uh(0,1:Ny) = uhrecvL
            uhsendL = uh(1,1:Ny)
            call MPI_SEND(uhsendL,Ny,MPI_REAL8,Left,myid + Nx0*Ny0,MPI_COMM_WORLD,ierr)
        end if
        
    end if
    
    ! y-direction
    if (mod(myidy,2) == 1) then ! odd
        
        if (Up >= 0) then
            uhsendU = uh(1:Nx,Ny)
            call MPI_SEND(uhsendU,Nx,MPI_REAL8,Up,myid,MPI_COMM_WORLD,ierr)
            call MPI_RECV(uhrecvU,Nx,MPI_REAL8,Up,Up + Nx0*Ny0,MPI_COMM_WORLD,status,ierr)
            uh(1:Nx,Ny1) = uhrecvU
        end if
        
        if (Down >= 0) then
            uhsendD = uh(1:Nx,1)
            call MPI_SEND(uhsendD,Nx,MPI_REAL8,Down,myid + Nx0*Ny0,MPI_COMM_WORLD,ierr)
            call MPI_RECV(uhrecvD,Nx,MPI_REAL8,Down,Down,MPI_COMM_WORLD,status,ierr)
            uh(1:Nx,0) = uhrecvD
        end if

    else ! even
        
        if (Up >= 0) then
            call MPI_RECV(uhrecvU,Nx,MPI_REAL8,Up,Up + Nx0*Ny0,MPI_COMM_WORLD,status,ierr)
            uh(1:Nx,Ny1) = uhrecvU
            uhsendU = uh(1:Nx,Ny)
            call MPI_SEND(uhsendU,Nx,MPI_REAL8,Up,myid,MPI_COMM_WORLD,ierr)
        end if
        
        if (Down >= 0) then
            call MPI_RECV(uhrecvD,Nx,MPI_REAL8,Down,Down,MPI_COMM_WORLD,status,ierr)
            uh(1:Nx,0) = uhrecvD
            uhsendD = uh(1:Nx,1)
            call MPI_SEND(uhsendD,Nx,MPI_REAL8,Down,myid + Nx0*Ny0,MPI_COMM_WORLD,ierr)
        end if
        
    end if
        
    end subroutine set_bc1
                
                
                
                
    