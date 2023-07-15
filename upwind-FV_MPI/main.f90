    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny
    parameter(N_process = 16)
    parameter(Nx0 = 1024, Ny0 = 1024)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,pi
    parameter(pi = 4*atan(1d0))
    integer bcR,bcL,bcU,bcD
    
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

    tend = 2*pi
    
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
    
    !print *,"process",myid1,"is alive,the index is",myidx,myidy
    
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
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    call system_clock(t2)
    
    if (myid1 == 1) then
        open(unit = 1,file = 'time.txt')
        write(1,*) t2 - t1
        print *,t2 - t1
    end if
    
    end subroutine Euler_Forward
    
    !**********************************************************************
    
    subroutine Lh
    
    use com
    
    call set_bc
    
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
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
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
                
                
                
                
    