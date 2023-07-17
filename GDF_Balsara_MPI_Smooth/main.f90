    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP
    parameter(N_process = 16)
    parameter(Nx0 = 256, Ny0 = 256, Lphi = 0, kk = 2, NumEq = 8, NumGLP = 3, RKorder = 4, flux_type = 2)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) pi
    parameter(dimPk = (kk + 2)*(kk + 3)/2)
    parameter(dimPk1 = (kk + 1)*(kk + 2)/2)
    parameter(Nphi = max(2*Lphi - 1,0))
    parameter(Nphi1 = Nphi + 1)
    parameter(gamma = 5d0/3d0)
    parameter(gamma1 = gamma - 1)
    parameter(pi = 4*atan(1d0))
    
    ! The numerical solution and mesh
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,tRK,t1,t2,alphax,alphay
    real(8) uh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),du(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uI(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),uII(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uh00(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uh0(Nx0,Ny0,0:Nphi,dimPk,NumEq)
    real(8) hx,hy,Xc(Nx),Yc(Ny),Xc0(Nx0),Yc0(Ny0),Phi(0:Nphi),hx1,hy1,hphi
    real(8) Bx(0:Nx,0:Ny1,0:Nphi,kk + 1),By(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) Bx0(0:Nx0,Ny0,0:Nphi,kk + 1),By0(Nx0,0:Ny0,0:Nphi,kk + 1)
    real(8) dBx(0:Nx,0:Ny1,0:Nphi,kk + 1),dBy(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxI(0:Nx,0:Ny1,0:Nphi,kk + 1),ByI(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxII(0:Nx,0:Ny1,0:Nphi,kk + 1),ByII(0:Nx1,0:Ny,0:Nphi,kk + 1)
    
    ! The basis
    real(8) lambda(NumGLP),weight(NumGLP),sink(0:Nphi,Lphi),cosk(0:Nphi,Lphi)
    real(8) phiG(NumGLP,NumGLP,dimPk),phixG(NumGLP,NumGLP,dimPk),phiyG(NumGLP,NumGLP,dimPk),mm(dimPk)
    real(8) phiGR(NumGLP,dimPk), phiGL(NumGLP,dimPk), phiGU(NumGLP,dimPk), phiGD(NumGLP,dimPk)
    real(8) phiRU(dimPk), phiLU(dimPk), phiRD(dimPk), phiLD(dimPk)
    real(8) EzG(NumGLP,kk + 1),EzxG(NumGLP,kk + 1),EzyG(NumGLP,kk + 1),mmE(kk + 1)
    real(8) EzR(kk + 1),EzL(kk + 1),EzU(kk + 1),EzD(kk + 1)
    
    ! The Lh
    real(8) uGint3D(NumGLP,NumGLP,0:Nphi,NumEq),uGint(NumGLP,NumGLP,NumEq)
    real(8) RHSC(NumGLP,NumGLP,0:Nphi,NumEq),RHSCopen,RG(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHS(NumGLP,NumGLP,0:Nphi,NumEq),Fzsin(Lphi),Fzcos(Lphi),Fzzsin(Lphi),Fzzcos(Lphi)
    real(8) FR1(NumEq),FL1(NumEq),UR1(NumEq),UL1(NumEq),Fhat1(NumEq),SR,SL
    real(8) URstar(NumEq),ULstar(NumEq),Ustar(NumEq),UUstar(NumEq),UDstar(NumEq)
    real(8) URU1(NumEq),ULU1(NumEq),URD1(NumEq),ULD1(NumEq)
    real(8) URstarstar(NumEq),ULstarstar(NumEq),Ezhat
    real(8) EzRL(0:Nx,Ny,0:Nphi,NumGLP), EzUD(Nx,0:Ny,0:Nphi,NumGLP), EzVertex(0:Nx,0:Ny,0:Nphi)
    real(8) URU(0:Nx1,0:Ny1,0:Nphi,NumEq),ULU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URD(0:Nx1,0:Ny1,0:Nphi,NumEq),ULD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) L2(NumEq)
    
    integer bcR,bcL,bcU,bcD,direction
    integer myid,myid1,the_id,the_id2
    integer myidx,myidy,the_idx,the_idy
    integer numprocs, namelen, rc,ierr,status(MPI_STATUS_SIZE),myid0
    character * (MPI_MAX_PROCESSOR_NAME) processor_name
    
    end module com
    
    !*****************************************************************************************************
    
    module init1
    
    use com
    
    contains
    
    function r2(x,y)
    real(8) r2,x,y
    r2 = x**2 + y**2
    end function r2
    
    function p(x,y)
    real(8) x,y,p
    p = 1 - r2(x,y)/(8.0d0*pi**2)*exp(1 - r2(x,y))
    end function p
    
    function rho(x,y)
    real(8) x,y,rho
    rho = 1
    end function rho
    
    function v1(x,y)
    real(8) v1,x,y
    v1 = 1 - 1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y)))
    end function v1
    
    function v2(x,y)
    real(8) v2,x,y
    v2 = 1 + 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y)))
    end function v2
    
    function v3(x,y)
    real(8) v3,x,y
    v3 = 0
    end function v3
    
    function B1(x,y)
    real(8) B1,x,y
    B1 = -1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y)))
    end function B1
    
    function B2(x,y)
    real(8) B2,x,y
    B2 = 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y)))
    end function B2
    
    function B3(x,y)
    real(8) B3,x,y
    B3 = 0
    end function B3

    function U1(x,y,z)
    real(8) U1,x,y,z
    U1 = rho(x,y)
    end function U1
    
    function U2(x,y,z)
    real(8) U2,x,y,z
    U2 = rho(x,y)*v1(x,y)
    end function U2
    
    function U3(x,y,z)
    real(8) U3,x,y,z
    U3 = rho(x,y)*v2(x,y)
    end function U3
    
    function U4(x,y,z)
    real(8) U4,x,y,z
    U4 = rho(x,y)*v3(x,y)
    end function U4
    
    function U5(x,y,z)
    real(8) U5,x,y,z
    U5 = p(x,y)/gamma1 + 0.5d0*rho(x,y)*(v1(x,y)**2 + v2(x,y)**2 + v3(x,y)**2) + 0.5d0*(B1(x,y)**2 + B2(x,y)**2 + B3(x,y)**2)
    end function U5
    
    function U6(x,y,z)
    real(8) U6,x,y,z
    U6 = B1(x,y)
    end function U6
    
    function U7(x,y,z)
    real(8) U7,x,y,z
    U7 = B2(x,y)
    end function U7
    
    function U8(x,y,z)
    real(8) U8,x,y,z
    U8 = B3(x,y)
    end function U8
    
    subroutine mesh
    
    use com
    
    xa = -10
    xb = 10
    ya = -10
    yb = 10

    bcL = 1
    bcR = 1
    bcU = 1
    bcD = 1

    tend = 20
    
    end subroutine mesh
    
    end module init1
    
    !*****************************************************************************************************
    
    module init2
    
    use com
    
    contains
    
    function p(x,y)
    real(8) x,y,p
    p = gamma
    end function p
    
    function rho(x,y)
    real(8) x,y,rho
    rho = gamma**2
    end function rho
    
    function v1(x,y)
    real(8) v1,x,y
    v1 = -sin(y)
    end function v1
    
    function v2(x,y)
    real(8) v2,x,y
    v2 = sin(x)
    end function v2
    
    function v3(x,y)
    real(8) v3,x,y
    v3 = 0
    end function v3
    
    function B1(x,y)
    real(8) B1,x,y
    B1 = -sin(y)
    end function B1
    
    function B2(x,y)
    real(8) B2,x,y
    B2 = sin(2*x)
    end function B2
    
    function B3(x,y)
    real(8) B3,x,y
    B3 = 0
    end function B3

    function U1(x,y,z)
    real(8) U1,x,y,z
    U1 = rho(x,y)
    end function U1
    
    function U2(x,y,z)
    real(8) U2,x,y,z
    U2 = rho(x,y)*v1(x,y)
    end function U2
    
    function U3(x,y,z)
    real(8) U3,x,y,z
    U3 = rho(x,y)*v2(x,y)
    end function U3
    
    function U4(x,y,z)
    real(8) U4,x,y,z
    U4 = rho(x,y)*v3(x,y)
    end function U4
    
    function U5(x,y,z)
    real(8) U5,x,y,z
    U5 = p(x,y)/gamma1 + 0.5d0*rho(x,y)*(v1(x,y)**2 + v2(x,y)**2 + v3(x,y)**2) + 0.5d0*(B1(x,y)**2 + B2(x,y)**2 + B3(x,y)**2)
    end function U5
    
    function U6(x,y,z)
    real(8) U6,x,y,z
    U6 = B1(x,y)
    end function U6
    
    function U7(x,y,z)
    real(8) U7,x,y,z
    U7 = B2(x,y)
    end function U7
    
    function U8(x,y,z)
    real(8) U8,x,y,z
    U8 = B3(x,y)
    end function U8
    
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

    tend = 0.5
    
    end subroutine mesh
    
    end module init2
    
    !*****************************************************************************************************
    
    program main
    
    use com
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    call CPU_TIME(t1)
    
    myid1 = myid + 1
    
    myidx = mod(myid1,Nx_process)
    if (myidx == 0) then
        myidx = Nx_process
    end if
    
    myidy = (myid1 - myidx)/Nx_process + 1
    
    !print *,"process",myid1,"is alive,the index is",myidx,myidy
    
    call init_data
    
    if (RKorder == 1) then
        call Euler_Forward
    else if (RKorder == 3) then
        call RK3
    else if (RKorder == 4) then
        call RK4
    end if
    
    call set_bc
    
    call save_solution
    
    call calculate_L2_Error
    
    call CPU_TIME(t2)
    
    if (myid1 == 1) then
        open(unit = 1,file = 'time.txt')
        write(1,*) t2 - t1
        print *,"Run time is",t2 - t1,"second"
        close(1)
    end if
    
    call MPI_FINALIZE(rc)
    
    end program main
    
    !*****************************************************************************************************
    
    subroutine init_data
    
    use com
    
    !use init1
    use init2
    
    call mesh
    
    hx = (xb - xa)/Nx0
    hy = (yb - ya)/Ny0
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    hphi = 2d0*pi/Nphi1
    
    call get_basis
    
    uh0 = 0
    uh = 0
    
    do i = 1,Nx0
        Xc0(i) = xa + (i - 0.5)*hx
    end do
    
    do j = 1,Ny0
        Yc0(j) = ya + (j - 0.5)*hy
    end do
    
    do k = 0,Nphi
        Phi(k) = k*hphi
    end do
    
    ! L2 Pro for Uh
    do i = 1,Nx0
        do j = 1,Ny0
            do k = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh0(i,j,k,d,1) = uh0(i,j,k,d,1) + 0.25*weight(i1)*weight(j1)*U1(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,2) = uh0(i,j,k,d,2) + 0.25*weight(i1)*weight(j1)*U2(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,3) = uh0(i,j,k,d,3) + 0.25*weight(i1)*weight(j1)*U3(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,4) = uh0(i,j,k,d,4) + 0.25*weight(i1)*weight(j1)*U4(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,5) = uh0(i,j,k,d,5) + 0.25*weight(i1)*weight(j1)*U5(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,6) = uh0(i,j,k,d,6) + 0.25*weight(i1)*weight(j1)*U6(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,7) = uh0(i,j,k,d,7) + 0.25*weight(i1)*weight(j1)*U7(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh0(i,j,k,d,8) = uh0(i,j,k,d,8) + 0.25*weight(i1)*weight(j1)*U8(Xc0(i) + hx1*lambda(i1),Yc0(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        uh0(:,:,:,d,:) = uh0(:,:,:,d,:)/mm(d)
    end do
    
    ! L2 Pro for Bx and By
    Bx0 = 0
    By0 = 0
    do i = 0,Nx0
        do j = 1,Ny0
            do k = 0,Nphi
                do d = 1,kk + 1
                    do j1 = 1,NumGLP
                        Bx0(i,j,k,d) = Bx0(i,j,k,d) + 0.5*weight(j1)*EzG(j1,d)*U6(xa + i*hx,Yc0(j) + hy1*lambda(j1),Phi(k))
                    end do
                end do
            end do
        end do
    end do
    
    do i = 1,Nx0
        do j = 0,Ny0
            do k = 0,Nphi
                do d = 1,kk + 1
                    do i1 = 1,NumGLP
                        By0(i,j,k,d) = By0(i,j,k,d) + 0.5*weight(i1)*EzG(i1,d)*U7(Xc0(i) + hx1*lambda(i1),ya + j*hy,Phi(k))
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,kk + 1
        Bx0(:,:,:,d) = Bx0(:,:,:,d)/mmE(d)
        By0(:,:,:,d) = By0(:,:,:,d)/mmE(d)
    end do
    
    uh(0:Nx1,0:Ny1,:,:,:) = uh0((myidx - 1)*Nx:myidx*Nx + 1,(myidy - 1)*Ny:myidy*Ny + 1,:,:,:)
    
    Bx(0:Nx,1:Ny,:,:) = Bx0((myidx - 1)*Nx:myidx*Nx,(myidy - 1)*Ny + 1:myidy*Ny,:,:)
    By(1:Nx,0:Ny,:,:) = By0((myidx - 1)*Nx + 1:myidx*Nx,(myidy - 1)*Ny:myidy*Ny,:,:)
    
    end subroutine init_data
    
    !*****************************************************************************************************
    
    subroutine save_solution
    
    use com
    
    !call div_free_Balsara
    !uh = 0
    !uh(1:Nx,1:Ny,:,1,6) = Ezvertex(1:Nx,1:Ny,:)
    !uh(1:Nx,1:Ny,:,1,7) = Ezvertex(1:Nx,1:Ny,:)
    !uh(1:Nx,1:Ny,:,1,6) = dBx(1:Nx,1:Ny,:,1)
    !uh(1:Nx,1:Ny,:,1,7) = dBy(1:Nx,1:Ny,:,1)
    !uh(1:Nx,1:Ny,:,1,6) = dBx(1:Nx,1:Ny,:,1)
    !uh(1:Nx,1:Ny,:,1,7) = dBy(1:Nx,1:Ny,:,1)
    !uh(1:Nx,1:Ny,:,1,6) = dBx(1:Nx,1:Ny,:,3)
    !uh(1:Nx,1:Ny,:,1,7) = dBy(1:Nx,1:Ny,:,3)
    uh0 = 0
    uh0((myidx - 1)*Nx + 1:myidx*Nx,(myidy - 1)*Ny + 1:myidy*Ny,0:Nphi,1:dimPk,1:NumEq) = uh(1:Nx,1:Ny,0:Nphi,1:dimPk,1:NumEq)
    
    do the_id = 2,N_process
        
        i = mod(the_id,Nx_process)
        if (i == 0) then
            i = Nx_process
        end if
        j = (the_id - i)/Nx_process + 1
        
        do d = 1,dimPk
            do n = 1,NumEq
                do k = 0,Nphi
                    
                    if (myid1 == the_id) then
                        call MPI_SEND(uh(1:Nx,1:Ny,k,d,n),Nx*Ny,MPI_REAL8,0,the_id,MPI_COMM_WORLD,ierr)
                    end if
                    
                    if (myid1 == 1) then
                        call MPI_RECV(uh0((i - 1)*Nx + 1:i*Nx,(j - 1)*Ny + 1:j*Ny,k,d,n),Nx*Ny,MPI_REAL8,the_id - 1,the_id,MPI_COMM_WORLD,status,ierr)
                    end if
                    
                end do
            end do
        end do
                    
        
    end do
            
    if (myid1 == 1) then
        
        open(unit = 1,file = 'Q1.txt')
        open(unit = 2,file = 'Q2.txt')
        open(unit = 3,file = 'Q3.txt')
        open(unit = 4,file = 'Q4.txt')
        open(unit = 5,file = 'Q5.txt')
        open(unit = 6,file = 'Q6.txt')
        open(unit = 7,file = 'Q7.txt')
        open(unit = 8,file = 'Q8.txt')
        
        open(unit = 9,file = 'Xc.txt')
        open(unit = 10,file = 'Yc.txt')
        
        !do d = 1,dimPk
        !    do k = 0,Nphi
        !        do j = 1,Ny0
        !            do i = 1,Nx0
        !                write(1,*) uh0(i,j,k,d,1)
        !                write(2,*) uh0(i,j,k,d,2)
        !                write(3,*) uh0(i,j,k,d,3)
        !                write(4,*) uh0(i,j,k,d,4)
        !                write(5,*) uh0(i,j,k,d,5)
        !                write(6,*) uh0(i,j,k,d,6)
        !                write(7,*) uh0(i,j,k,d,7)
        !                write(8,*) uh0(i,j,k,d,8)
        !            end do
        !        end do
        !    end do
        !end do
        
        do d = 1,dimPk
            do j = 1,Ny0
                do i = 1,Nx0
                    write(1,*) uh0(i,j,0,d,1)
                    write(2,*) uh0(i,j,0,d,2)
                    write(3,*) uh0(i,j,0,d,3)
                    write(4,*) uh0(i,j,0,d,4)
                    write(5,*) uh0(i,j,0,d,5)
                    write(6,*) uh0(i,j,0,d,6)
                    write(7,*) uh0(i,j,0,d,7)
                    write(8,*) uh0(i,j,0,d,8)
                end do
            end do
        end do
        
        do i = 1,Nx0
            write(9,*) Xc0(i)
        end do
        
        do j = 1,Ny0
            write(10,*) Yc0(j)
        end do
        
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
        close(9)
        close(10)
        
    end if
    
    end subroutine save_solution
    
    !*****************************************************************************************************
    
    subroutine get_basis
    
    use com
    
    ! get Gauss points
    if (NumGLP == 2) then
        lambda(1) = -0.5773502691896257645091488
        lambda(2) = 0.5773502691896257645091488
        
        weight(1) = 1
        weight(2) = 1
    else if (NumGLP == 3) then
        lambda(1) = -0.7745966692414833770358531d0
        lambda(2) = 0
        lambda(3) = 0.7745966692414833770358531d0
        
        weight(1) = 0.5555555555555555555555556d0
        weight(2) = 0.8888888888888888888888889d0
        weight(3) = 0.5555555555555555555555556d0
    else if (NumGLP == 4) then
        lambda(1) = -0.8611363115940525752239465d0
        lambda(2) = -0.3399810435848562648026658d0
        lambda(3) = 0.3399810435848562648026658d0
        !lambda(4) = 0.8611363115940525752239465d0
        
        weight(1) = 0.3478548451374538573730639d0
        weight(2) = 0.6521451548625461426269361d0
        weight(3) = 0.6521451548625461426269361d0   
        !weight(4) = 0.3478548451374538573730639d0
    else if (NumGLP == 5) then
        lambda(1) = -0.9061798459386639927976269     
        lambda(2) = -0.5384693101056830910363144     
        lambda(3) = 0                                 
        !lambda(4) = 0.5384693101056830910363144     
        !lambda(5) = 0.9061798459386639927976269     
        
        weight(1) = 0.2369268850561890875142640
        weight(2) = 0.4786286704993664680412915
        weight(3) = 0.5688888888888888888888889
        !weight(4) = 0.4786286704993664680412915
        !weight(5) = 0.2369268850561890875142640
    else if (NumGLP == 6) then
        lambda(1) = -0.9324695142031520278123016     
        lambda(2) = -0.6612093864662645136613996    
        !lambda(3) = -0.2386191860831969086305017     
        !lambda(4) = 0.2386191860831969086305017     
        !lambda(5) = 0.6612093864662645136613996     
        !lambda(6) = 0.9324695142031520278123016     
        
        weight(1) = 0.1713244923791703450402961
        weight(2) = 0.3607615730481386075698335
        !weight(3) = 0.4679139345726910473898703
        !weight(4) = 0.4679139345726910473898703
        !weight(5) = 0.3607615730481386075698335
        !weight(6) = 0.1713244923791703450402961
    end if
    
    phiRU(1) = 1
    phiRU(2) = 1
    phiRU(3) = 1
    phiRU(4) = 2d0/3d0
    phiRU(5) = 1
    phiRU(6) = 2d0/3d0
    phiRU(7) = 2d0/5d0
    phiRU(8) = 2d0/3d0
    phiRU(9) = 2d0/3d0
    phiRU(10) = 2d0/5d0
    
    phiLU(1) = 1
    phiLU(2) = -1
    phiLU(3) = 1
    phiLU(4) = 2d0/3d0
    phiLU(5) = -1
    phiLU(6) = 2d0/3d0
    phiLU(7) = -2d0/5d0
    phiLU(8) = 2d0/3d0
    phiLU(9) = -2d0/3d0
    phiLU(10) = 2d0/5d0
    
    phiRD(1) = 1
    phiRD(2) = 1
    phiRD(3) = -1
    phiRD(4) = 2d0/3d0
    phiRD(5) = -1
    phiRD(6) = 2d0/3d0
    phiRD(7) = 2d0/5d0
    phiRD(8) = -2d0/3d0
    phiRD(9) = 2d0/3d0
    phiRD(10) = -2d0/5d0
    
    phiLD(1) = 1
    phiLD(2) = -1
    phiLD(3) = -1
    phiLD(4) = 2d0/3d0
    phiLD(5) = 1
    phiLD(6) = 2d0/3d0
    phiLD(7) = -2d0/5d0
    phiLD(8) = -2d0/3d0
    phiLD(9) = -2d0/3d0
    phiLD(10) = -2d0/5d0
    
    do i = 1,NumGLP
        do j = 1,NumGLP
            phiG(i,j,1) = 1
            phiGR(j,1) = 1
            phiGL(j,1) = 1
            phiGU(i,1) = 1
            phiGD(i,1) = 1
            phixG(i,j,1) = 0
            phiyG(i,j,1) = 0
            mm(1) = 1
            
            phiG(i,j,2) = lambda(i)
            phiGR(j,2) = 1
            phiGL(j,2) = -1
            phiGU(i,2) = lambda(i)
            phiGD(i,2) = lambda(i)
            phixG(i,j,2) = 1d0/hx1
            phiyG(i,j,2) = 0
            mm(2) = 1d0/3d0
                
            phiG(i,j,3) = lambda(j)
            phiGR(j,3) = lambda(j)
            phiGL(j,3) = lambda(j)
            phiGU(i,3) = 1
            phiGD(i,3) = -1
            phixG(i,j,3) = 0
            phiyG(i,j,3) = 1d0/hy1
            mm(3) = 1d0/3d0
                
            phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
            phiGR(j,4) = 2d0/3d0
            phiGL(j,4) = 2d0/3d0
            phiGU(i,4) = lambda(i)**2 - 1d0/3d0
            phiGD(i,4) = lambda(i)**2 - 1d0/3d0
            phixG(i,j,4) = 2d0*lambda(i)/hx1
            phiyG(i,j,4) = 0
            mm(4) = 4d0/45d0
                    
            phiG(i,j,5) = lambda(i)*lambda(j)
            phiGR(j,5) = lambda(j)
            phiGL(j,5) = -lambda(j)
            phiGU(i,5) = lambda(i)
            phiGD(i,5) = -lambda(i)
            phixG(i,j,5) = lambda(j)/hx1
            phiyG(i,j,5) = lambda(i)/hy1
            mm(5) = 1d0/9d0
                    
            phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
            phiGR(j,6) = lambda(j)**2 - 1d0/3d0
            phiGL(j,6) = lambda(j)**2 - 1d0/3d0
            phiGU(i,6) = 2d0/3d0
            phiGD(i,6) = 2d0/3d0
            phixG(i,j,6) = 0
            phiyG(i,j,6) = 2d0*lambda(j)/hy1
            mm(6) = 4d0/45d0
                
            phiG(i,j,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGR(j,7) = 2d0/5d0
            phiGL(j,7) = -2d0/5d0
            phiGU(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGD(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phixG(i,j,7) = (3*lambda(i)**2 - 3d0/5d0)/hx1
            phiyG(i,j,7) = 0
            mm(7) = 4d0/175d0
                
            phiG(i,j,8) = (lambda(i)**2 - 1d0/3d0)*(lambda(j))
            phiGR(j,8) = (2d0/3d0)*(lambda(j))
            phiGL(j,8) = (2d0/3d0)*(lambda(j))
            phiGU(i,8) = (lambda(i)**2 - 1d0/3d0)
            phiGD(i,8) = -(lambda(i)**2 - 1d0/3d0)
            phixG(i,j,8) = 2d0*lambda(i)*lambda(j)/hx1
            phiyG(i,j,8) = (lambda(i)**2 - 1d0/3d0)/hy1
            mm(8) = 4d0/135d0
                
            phiG(i,j,9) = (lambda(i))*(lambda(j)**2 - 1d0/3d0)
            phiGR(j,9) = (lambda(j)**2 - 1d0/3d0)
            phiGL(j,9) = -(lambda(j)**2 - 1d0/3d0)
            phiGU(i,9) = lambda(i)*(2d0/3d0)
            phiGD(i,9) = lambda(i)*(2d0/3d0)
            phixG(i,j,9) = (lambda(j)**2 - 1d0/3d0)/hx1
            phiyG(i,j,9) = 2d0*lambda(i)*lambda(j)/hy1
            mm(9) = 4d0/135d0
                
            phiG(i,j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGR(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGL(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGU(i,10) = 2d0/5d0
            phiGD(i,10) = -2d0/5d0
            phixG(i,j,10) = 0
            phiyG(i,j,10) = (3*lambda(j)**2 - 3d0/5d0)/hy1
            mm(10) = 4d0/175d0
        end do
        
        EzG(i,1) = 1
        EzG(i,2) = lambda(i)
        EzG(i,3) = lambda(i)**2 - 1d0/3d0
        
        EzxG(i,1) = 0
        EzxG(i,2) = 1/hx1
        EzxG(i,3) = 2*lambda(i)/hx1
        
        EzyG(i,1) = 0
        EzyG(i,2) = 1/hy1
        EzyG(i,3) = 2*lambda(i)/hy1
        
    end do
    
    EzR(1) = 1
    EzR(2) = 1
    EzR(3) = 2d0/3d0
    
    EzL(1) = 1
    EzL(2) = -1
    EzL(3) = 2d0/3d0
    
    EzU = EzR
    EzD = EzL
    
    mmE(1) = 1
    mmE(2) = 1d0/3d0
    mmE(3) = 4d0/45d0
    
    do k = 1,Lphi
        do i = 0,Nphi
            sink(i,k) = sin(k*Phi(i))
            cosk(i,k) = cos(k*Phi(i))
        end do
    end do
    
    end subroutine get_basis
    
    !*****************************************************************************************************
    
    subroutine calculate_L2_Error
    
    use com
    
    use init1
    
    if (tend > 3) then
        tend = 0
    end if
    
    L2 = 0
    
    if (myid1 == 1) then
        
        do i = 1,Nx0
            do j = 1,Ny0
                do k = 0,Nphi
                
                    uGint = 0
                    do d = 1,dimPk
                        do n = 1,NumEq
                            uGint(:,:,n) = uGint(:,:,n) + uh0(i,j,k,d,n)*phiG(:,:,d)
                        end do
                    end do
                
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            L2(1) = L2(1) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,1) - U1(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(2) = L2(2) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,2) - U2(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(3) = L2(3) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,3) - U3(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(4) = L2(4) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,4) - U4(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(5) = L2(5) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,5) - U5(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(6) = L2(6) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,6) - U6(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(7) = L2(7) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,7) - U7(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                            L2(8) = L2(8) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,8) - U8(Xc0(i) + hx1*lambda(i1) - tend,Yc0(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        end do
                    end do
                
                end do
            end do
        end do
        
    end if
    
    L2 = (L2/(Nx0*Ny0*Nphi1))**0.5d0
    
    if (myid1 == 1) then
        print *,"The L2 Error:"
        print *,"rho    :",L2(1)
        print *,"rho u1 :",L2(2)
        print *,"rho u2 :",L2(3)
        print *,"rho u3 :",L2(4)
        print *,"E      :",L2(5)
        print *,"B1     :",L2(6)
        print *,"B2     :",L2(7)
        print *,"B3     :",L2(8)
    end if
    
    end subroutine calculate_L2_Error
    
    !*****************************************************************************************************
    
    subroutine Euler_Forward
    
    use com
    
    CFL = 0.02
    t = 0
    
    call div_free_Balsara
    
    do while (t < tend)
        
        call calculate_dt
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        call Lh
        
        uh = uh + dt*du
        Bx = Bx + dt*dBx
        By = By + dt*dBy
        
        call div_free_Balsara
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    end subroutine Euler_Forward
    
    !*****************************************************************************************************
    
    subroutine RK3
    
    use com
    
    CFL = 0.12
    t = 0
    
    do while (t < tend)
        
        call calculate_dt
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        ! Stage I
        call Lh
        
        uh00 = uh
        
        uI = uh + dt*du
        
        uh = uI
        
        ! Stage II
        call Lh
        
        uII = (3d0/4d0)*uh00 + (1d0/4d0)*uh + (1d0/4d0)*dt*du
        
        uh = uII
        
        ! Stage III
        call Lh
        
        uh = (1d0/3d0)*uh00 + (2d0/3d0)*uh + (2d0/3d0)*dt*du
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    end subroutine RK3
    
    !*****************************************************************************************************
    
    subroutine RK4
    
    use com
    
    CFL = 1.4
    t = 0
    
    call div_free_Balsara
    
    do while (t < tend)
        
        call calculate_dt
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        uI = uh
        BxI = Bx
        ByI = By
        
        uII = uh
        BxII = Bx
        ByII = By
        
        ! Stage I
        do i = 1,5
            
            call Lh
            
            uI = uh + (dt/6d0)*du
            BxI = Bx + (dt/6d0)*dBx
            ByI = By + (dt/6d0)*dBy
            
            uh = uI
            Bx = BxI
            By = ByI
            
            call div_free_Balsara
            
        end do
        
        uII = 0.04d0*uII + 0.36d0*uI
        BxII = 0.04d0*BxII + 0.36d0*BxI
        ByII = 0.04d0*ByII + 0.36d0*ByI
        
        uI = 15*uII - 5*uI
        BxI = 15*BxII - 5*BxI
        ByI = 15*ByII - 5*ByI
        
        uh = uI
        Bx = BxI
        By = ByI
        
        ! Stage II
        do i = 6,9
            
            call Lh
            
            uI = uh + (dt/6d0)*du
            BxI = Bx + (dt/6d0)*dBx
            ByI = By + (dt/6d0)*dBy
            
            uh = uI
            Bx = BxI
            By = ByI
            
            call div_free_Balsara
            
        end do
        
        call Lh
        
        uh = uII + 0.6d0*uI + (dt/10d0)*du
        Bx = BxII + 0.6d0*BxI + (dt/10d0)*dBx
        By = ByII + 0.6d0*ByI + (dt/10d0)*dBy
        
        call div_free_Balsara
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    end subroutine RK4
    
    !*****************************************************************************************************
    
    subroutine calculate_umax
    
    use com
    
    umax = 0
    umax1 = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                if (abs(uh(i,j,k,1,1)) > umax) then
                    umax = abs(uh(i,j,k,1,1))
                end if
            end do
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
    
    !*****************************************************************************************************
    
    subroutine calculate_dt
    
    use com
    
    alphax = 0
    alphay = 0
    alphax0 = 0
    alphay0 = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                call eigenvalueMm(alpha1,alpha2,uh(i,j,k,1,1),uh(i,j,k,1,2),uh(i,j,k,1,3),uh(i,j,k,1,4),uh(i,j,k,1,5),uh(i,j,k,1,6),uh(i,j,k,1,7),uh(i,j,k,1,8),1,0)
                if (abs(alpha1) > alphax .or. abs(alpha2) > alphax) then
                    alphax = max(abs(alpha1),abs(alpha2))
                end if
                call eigenvalueMm(alpha1,alpha2,uh(i,j,k,1,1),uh(i,j,k,1,2),uh(i,j,k,1,3),uh(i,j,k,1,4),uh(i,j,k,1,5),uh(i,j,k,1,6),uh(i,j,k,1,7),uh(i,j,k,1,8),0,1)
                if (abs(alpha1) > alphay .or. abs(alpha2) > alphay) then
                    alphay = max(abs(alpha1),abs(alpha2))
                end if
            end do
        end do
    end do
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    do the_id = 2,N_process
        
        if (myid1 == the_id) then
            call MPI_SEND(alphax,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(alphax0,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            if (alphax0 > alphax) then
                alphax = alphax0
            end if
        end if
        
        if (myid1 == the_id) then
            call MPI_SEND(alphay,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(alphay0,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            if (alphay0 > alphay) then
                alphay = alphay0
            end if
        end if
        
    end do
    
    do the_id = 2,N_process
        
        if (myid1 == 1) then
            call MPI_SEND(alphax,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == the_id) then
            call MPI_RECV(alphax,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_SEND(alphay,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == the_id) then
            call MPI_RECV(alphay,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
        end if
        
    end do
    
    dt = CFL/(alphax/hx + alphay/hy)
    
    end subroutine calculate_dt

    !*****************************************************************************************************
    
    subroutine eigenvalueMm(Amax,Amin,rho,rhou,rhov,rhow,E,B1,B2,B3,n1,n2)
    
    use com
    
    real(8) u,v,w,p,c,BP,Bn,un,cf,n3
    
    n3 = 0
    
    u = rhou/rho
    v = rhov/rho
    w = rhow/rho
    
    BP = B1**2 + B2**2 + B3**2
    Bn = B1*n1 + B2*n2 + B3*n3
    un = u*n1 + v*n2 + w*n3
    
    p = gamma1*(E - 0.5d0*rho*(u**2 + v**2 + w**2) - 0.5d0*BP)
    
    c = sqrt(abs(gamma*p/rho))
    
    cf = sqrt(abs( 0.5d0*(c**2 + BP/rho + sqrt((c**2 + BP/rho)**2 - 4*c**2*Bn**2/rho) ) ))
    
    Amax = un + cf
    Amin = un - cf
    
    end subroutine eigenvalueMm
    
    !*****************************************************************************************************
    
    subroutine set_bc
    
    use com
    
    do i = 1,Nx_process
        do j = 1,Ny_process
            
            the_id = i + Nx_process*(j - 1)
            
            do k = 0,Nphi
                do n = 1,NumEq
                    do d = 1,dimPk
                        
                        ! The Right condition
                        if (i == Nx_process) then
                            the_idx = 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcR == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                                end if
                            end if
            
                        else
                
                            the_idx = i + 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                        end if
                        
                        ! The Left condition
                        if (i == 1) then
                            the_idx = Nx_process
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcL == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i - 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                        end if
            
                        ! The Up condition
                        if (j == Ny_process) then
                            the_idx = i
                            the_idy = 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcU == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j + 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,Ny1,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                        end if
            
                        ! The Down condition
                        if (j == 1) then
                            the_idx = i
                            the_idy = Ny_process
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcD == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j - 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,0,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                        end if
                        
                    end do
                end do
            end do
            
        end do
    end do
    
    end subroutine set_bc
    
    !*****************************************************************************************************
    
    subroutine Lh
    
    use com
    
    real(8) Fx(NumGLP,NumGLP,0:Nphi,NumEq), Fy(NumGLP,NumGLP,0:Nphi,NumEq), Fz(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) rhoij,uij,vij,wij,Eij,B1ij,B2ij,B3ij,pij,Sij,Tij,Kij,rB1ij,rB2ij,rB3ij
    
    real(8),allocatable :: UR(:,:,:,:,:),UL(:,:,:,:,:),UU(:,:,:,:,:),UD(:,:,:,:,:)
    real(8),allocatable :: FR(:,:,:,:,:),FL(:,:,:,:,:),FU(:,:,:,:,:),FD(:,:,:,:,:)
    !real(8),allocatable :: URU(:,:,:,:),ULU(:,:,:,:),URD(:,:,:,:),ULD(:,:,:,:)
    !real(8),allocatable :: EzVertex(:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:)
    !real(8),allocatable :: EzR1(:,:,:,:),EzL1(:,:,:,:),ErU(:,:,:,:),ErD(:,:,:,:)
    !real(8),allocatable :: EphiRL(:,:,:,:), EzRL(:,:,:,:), EphiUD(:,:,:,:), ErUD(:,:,:,:)
    real(8),allocatable :: RHS1(:,:,:,:), RHS2(:,:,:,:)
    
    allocate(UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(UD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(FR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(FD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    !allocate(URU(0:Nx1,0:Ny1, 0:Nphi, NumEq))
    !allocate(ULU(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(URD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(ULD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(EzVertex(0:Nx,0:Ny,0:Nphi))
    
    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    
    !allocate(EzR1(0:Nx,Ny,0:Nphi,NumGLP))
    !allocate(EzL1(1:Nx1,Ny,0:Nphi,NumGLP))
    !allocate(ErU(Nx,0:Ny,0:Nphi,NumGLP))
    !allocate(ErD(Nx,0:Ny,0:Nphi,NumGLP))
    
    !allocate(EphiRL(0:Nx,Ny,0:Nphi,NumGLP))
    !allocate(EzRL(0:Nx,Ny,0:Nphi,NumGLP))
    !allocate(EphiUD(Nx,0:Ny,0:Nphi,NumGLP))
    !allocate(ErUD(Nx,0:Ny,0:Nphi,NumGLP))
    
    allocate(RHS1(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(RHS2(Nx,0:Ny,0:Nphi,NumGLP))
    
    call set_bc
    call set_bc
    
    du = 0
    
    ! calculate the Volume integral
    do i = 1,Nx
        do j = 1,Ny
            
            uGint3D = 0
            RHS = 0
            do k = 0,Nphi
                do n = 1,NumEq
                    do d = 1,dimPk
                        uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do k = 0,Nphi
                        rhoij = uGint3D(i1,j1,k,1)
                        uij = uGint3D(i1,j1,k,2)/rhoij
                        vij = uGint3D(i1,j1,k,3)/rhoij
                        wij = uGint3D(i1,j1,k,4)/rhoij
                        Eij = uGint3D(i1,j1,k,5)
                        B1ij = uGint3D(i1,j1,k,6)
                        B2ij = uGint3D(i1,j1,k,7)
                        B3ij = uGint3D(i1,j1,k,8)
    
                        pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2 + wij**2) - 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2))
    
                        Sij = pij + 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2)
                        Tij = Eij + Sij
                        Kij = uij*B1ij + vij*B2ij + wij*B3ij
    
                        Fx(i1,j1,k,1) = rhoij*uij
                        Fx(i1,j1,k,2) = rhoij*uij**2 + Sij - B1ij**2
                        Fx(i1,j1,k,3) = rhoij*uij*vij - B1ij*B2ij
                        Fx(i1,j1,k,4) = rhoij*uij*wij - B1ij*B3ij
                        Fx(i1,j1,k,5) = Tij*uij - Kij*B1ij
                        Fx(i1,j1,k,6) = 0
                        Fx(i1,j1,k,7) = uij*B2ij - vij*B1ij
                        Fx(i1,j1,k,8) = uij*B3ij - wij*B1ij
    
                        Fy(i1,j1,k,1) = rhoij*vij
                        Fy(i1,j1,k,2) = rhoij*uij*vij - B1ij*B2ij
                        Fy(i1,j1,k,3) = rhoij*vij**2 + Sij - B2ij**2
                        Fy(i1,j1,k,4) = rhoij*vij*wij - B2ij*B3ij
                        Fy(i1,j1,k,5) = Tij*vij - Kij*B2ij
                        Fy(i1,j1,k,6) = vij*B1ij - uij*B2ij
                        Fy(i1,j1,k,7) = 0
                        Fy(i1,j1,k,8) = vij*B3ij - wij*B2ij
                
                        Fz(i1,j1,k,1) = rhoij*wij
                        Fz(i1,j1,k,2) = rhoij*uij*wij - B1ij*B3ij
                        Fz(i1,j1,k,3) = rhoij*vij*wij - B2ij*B3ij
                        Fz(i1,j1,k,4) = rhoij*wij**2 + Sij - B3ij**2
                        Fz(i1,j1,k,5) = Tij*wij - Kij*B3ij
                        Fz(i1,j1,k,6) = wij*B1ij - uij*B3ij
                        Fz(i1,j1,k,7) = wij*B2ij - vij*B3ij
                        Fz(i1,j1,k,8) = 0
                    end do
                end do
            end do
            
            do d = 1,dimPk1
                do n = 1,NumEq
                    do k = 0,Nphi
                        do i1 = 1,NumGLP
                            do j1 = 1,NumGLP
                                if (d > 1) then
                                    du(i,j,k,d,n) = du(i,j,k,d,n) + 0.25d0*weight(i1)*weight(j1)*(Fx(i1,j1,k,n)*phixG(i1,j1,d) + Fy(i1,j1,k,n)*phiyG(i1,j1,d))
                                end if
                                !du(i,j,kk,d,n) = du(i,j,kk,d,n) + 0.25d0*weight(i1)*weight(j1)*RHS(i1,j1,kk,n)*phiG(i1,j1,d)
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    
    ! calculate the Numerical flux
    
    ! The x-flux
    UR = 0
    UL = 0
    
    do i = 0,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk
                    do n = 1,NumEq
                        UR(i,j,k,:,n) = UR(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR(:,d)
                        UL(i + 1,j,k,:,n) = UL(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do j1 = 1,NumGLP
                    rhoij = uR(i,j,k,j1,1)
                    uij = uR(i,j,k,j1,2)/rhoij
                    vij = uR(i,j,k,j1,3)/rhoij
                    wij = uR(i,j,k,j1,4)/rhoij
                    Eij = uR(i,j,k,j1,5)
                    B1ij = uR(i,j,k,j1,6)
                    B2ij = uR(i,j,k,j1,7)
                    B3ij = uR(i,j,k,j1,8)
    
                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2 + wij**2) - 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2))
    
                    Sij = pij + 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2)
                    Tij = Eij + Sij
                    Kij = uij*B1ij + vij*B2ij + wij*B3ij
    
                    FR(i,j,k,j1,1) = rhoij*uij
                    FR(i,j,k,j1,2) = rhoij*uij**2 + Sij - B1ij**2
                    FR(i,j,k,j1,3) = rhoij*uij*vij - B1ij*B2ij
                    FR(i,j,k,j1,4) = rhoij*uij*wij - B1ij*B3ij
                    FR(i,j,k,j1,5) = Tij*uij - Kij*B1ij
                    FR(i,j,k,j1,6) = 0
                    FR(i,j,k,j1,7) = uij*B2ij - vij*B1ij
                    FR(i,j,k,j1,8) = uij*B3ij - wij*B1ij
                end do
            end do
        end do
    end do
    
    do i = 1,Nx1
        do j = 1,Ny
            do k = 0,Nphi
                do j1 = 1,NumGLP
                    rhoij = uL(i,j,k,j1,1)
                    uij = uL(i,j,k,j1,2)/rhoij
                    vij = uL(i,j,k,j1,3)/rhoij
                    wij = uL(i,j,k,j1,4)/rhoij
                    Eij = uL(i,j,k,j1,5)
                    B1ij = uL(i,j,k,j1,6)
                    B2ij = uL(i,j,k,j1,7)
                    B3ij = uL(i,j,k,j1,8)
    
                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2 + wij**2) - 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2))
    
                    Sij = pij + 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2)
                    Tij = Eij + Sij
                    Kij = uij*B1ij + vij*B2ij + wij*B3ij
    
                    FL(i,j,k,j1,1) = rhoij*uij
                    FL(i,j,k,j1,2) = rhoij*uij**2 + Sij - B1ij**2
                    FL(i,j,k,j1,3) = rhoij*uij*vij - B1ij*B2ij
                    FL(i,j,k,j1,4) = rhoij*uij*wij - B1ij*B3ij
                    FL(i,j,k,j1,5) = Tij*uij - Kij*B1ij
                    FL(i,j,k,j1,6) = 0
                    FL(i,j,k,j1,7) = uij*B2ij - vij*B1ij
                    FL(i,j,k,j1,8) = uij*B3ij - wij*B1ij
                end do
            end do
        end do
    end do
    
    ! The y-Flux
    UU = 0
    UD = 0
    
    do i = 1,Nx
        do j = 0,Ny
            do k = 0,Nphi
                do d = 1,dimPk
                    do n = 1,NumEq
                        UU(i,j,k,:,n) = UU(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU(:,d)
                        UD(i,j + 1,k,:,n) = UD(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do k = 0,Nphi
                do i1 = 1,NumGLP
                    rhoij = UU(i,j,k,i1,1)
                    uij = UU(i,j,k,i1,2)/rhoij
                    vij = UU(i,j,k,i1,3)/rhoij
                    wij = UU(i,j,k,i1,4)/rhoij
                    Eij = UU(i,j,k,i1,5)
                    B1ij = UU(i,j,k,i1,6)
                    B2ij = UU(i,j,k,i1,7)
                    B3ij = UU(i,j,k,i1,8)
    
                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2 + wij**2) - 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2))
    
                    Sij = pij + 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2)
                    Tij = Eij + Sij
                    Kij = uij*B1ij + vij*B2ij + wij*B3ij
    
                    FU(i,j,k,i1,1) = rhoij*vij
                    FU(i,j,k,i1,2) = rhoij*uij*vij - B1ij*B2ij
                    FU(i,j,k,i1,3) = rhoij*vij**2 + Sij - B2ij**2
                    FU(i,j,k,i1,4) = rhoij*vij*wij - B2ij*B3ij
                    FU(i,j,k,i1,5) = Tij*vij - Kij*B2ij
                    FU(i,j,k,i1,6) = vij*B1ij - uij*B2ij
                    FU(i,j,k,i1,7) = 0
                    FU(i,j,k,i1,8) = vij*B3ij - wij*B2ij
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny1
            do k = 0,Nphi
                do i1 = 1,NumGLP
                    rhoij = UD(i,j,k,i1,1)
                    uij = UD(i,j,k,i1,2)/rhoij
                    vij = UD(i,j,k,i1,3)/rhoij
                    wij = UD(i,j,k,i1,4)/rhoij
                    Eij = UD(i,j,k,i1,5)
                    B1ij = UD(i,j,k,i1,6)
                    B2ij = UD(i,j,k,i1,7)
                    B3ij = UD(i,j,k,i1,8)
    
                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2 + wij**2) - 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2))
    
                    Sij = pij + 0.5d0*(B1ij**2 + B2ij**2 + B3ij**2)
                    Tij = Eij + Sij
                    Kij = uij*B1ij + vij*B2ij + wij*B3ij
    
                    FD(i,j,k,i1,1) = rhoij*vij
                    FD(i,j,k,i1,2) = rhoij*uij*vij - B1ij*B2ij
                    FD(i,j,k,i1,3) = rhoij*vij**2 + Sij - B2ij**2
                    FD(i,j,k,i1,4) = rhoij*vij*wij - B2ij*B3ij
                    FD(i,j,k,i1,5) = Tij*vij - Kij*B2ij
                    FD(i,j,k,i1,6) = vij*B1ij - uij*B2ij
                    FD(i,j,k,i1,7) = 0
                    FD(i,j,k,i1,8) = vij*B3ij - wij*B2ij
                end do
            end do
        end do
    end do
    
    ! calculate Fx hat
    do i = 0,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do j1 = 1,NumGLP
                    call eigenvalueMm(SRmax,SRmin,UR(i,j,k,j1,1),UR(i,j,k,j1,2),UR(i,j,k,j1,3),UR(i,j,k,j1,4),UR(i,j,k,j1,5),UR(i,j,k,j1,6),UR(i,j,k,j1,7),UR(i,j,k,j1,8),1,0)
                    call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,k,j1,1),UL(i + 1,j,k,j1,2),UL(i + 1,j,k,j1,3),UL(i + 1,j,k,j1,4),UL(i + 1,j,k,j1,5),UL(i + 1,j,k,j1,6),UL(i + 1,j,k,j1,7),UL(i + 1,j,k,j1,8),1,0)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FL(i + 1,j,k,j1,:)
                    FL1 = FR(i,j,k,j1,:)
                    UR1 = UL(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 1
                        call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 1
                        call HLLD_Flux
                    end if
                    Fxhat(i,j,k,j1,:) = Fhat1
                end do
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i = 1,Nx
        do j = 0,Ny
            do k = 0,Nphi
                do i1 = 1,NumGLP
                    call eigenvalueMm(SRmax,SRmin,UU(i,j,k,i1,1),UU(i,j,k,i1,2),UU(i,j,k,i1,3),UU(i,j,k,i1,4),UU(i,j,k,i1,5),UU(i,j,k,i1,6),UU(i,j,k,i1,7),UU(i,j,k,i1,8),0,1)
                    call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,k,i1,1),UD(i,j + 1,k,i1,2),UD(i,j + 1,k,i1,3),UD(i,j + 1,k,i1,4),UD(i,j + 1,k,i1,5),UD(i,j + 1,k,i1,6),UD(i,j + 1,k,i1,7),UD(i,j + 1,k,i1,8),0,1)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FD(i,j + 1,k,i1,:)
                    FL1 = FU(i,j,k,i1,:)
                    UR1 = UD(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 2
                        call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 2
                        call HLLD_Flux
                    end if
                    Fyhat(i,j,k,i1,:) = Fhat1
                end do
            end do
        end do
    end do
    
    ! calculate the Surface integral
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        do j1 = 1,NumGLP
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hx)*weight(j1)*(Fxhat(i,j,k,j1,n)*phiGR(j1,d) - Fxhat(i - 1,j,k,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        do i1 = 1,NumGLP
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hy)*weight(i1)*(Fyhat(i,j,k,i1,n)*phiGU(i1,d) - Fyhat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        du(:,:,:,d,:) = du(:,:,:,d,:)/mm(d)
    end do
    
    
    
    ! Scheme for Bx and By on the edge
    EzRL = -Fxhat(:,:,:,:,7)
    EzUD = Fyhat(:,:,:,:,6)
    
    ! calculate each Uh at vertex
    URU = 0
    ULU = 0
    URD = 0
    ULD = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do k = 0,Nphi
                do d = 1,dimPk
                    URU(i,j,k,:) = URU(i,j,k,:) + uh(i,j,k,d,:)*phiRU(d)
                    ULU(i,j,k,:) = ULU(i,j,k,:) + uh(i,j,k,d,:)*phiLU(d)
                    URD(i,j,k,:) = URD(i,j,k,:) + uh(i,j,k,d,:)*phiRD(d)
                    ULD(i,j,k,:) = ULD(i,j,k,:) + uh(i,j,k,d,:)*phiLD(d)
                end do
            end do
        end do
    end do
    
    ! calculate Ez at vertex
    do i = 0,Nx
        do j = 0,Ny
            do k = 0,Nphi
                URU1 = ULD(i + 1,j + 1,k,:)
                ULU1 = URD(i,j + 1,k,:)
                URD1 = ULU(i + 1,j,k,:)
                ULD1 = URU(i,j,k,:)
            
                if (flux_type == 1) then
                    call LF_Flux_2D
                else if (flux_type == 2) then
                    call HLL_Flux_2D
                else if (flux_type == 3) then
                    call HLLC_Flux_2D
                else if (flux_type == 4) then
                    call HLLD_Flux_2D
                end if
            
                EzVertex(i,j,k) = Ezhat
            end do
        end do
    end do
    
    dBx = 0
    dBy = 0
    
    ! DG scheme of Ez1
    do i = 0,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,kk + 1
                    do j1 = 1,NumGLP
                        dBx(i,j,k,d) = dBx(i,j,k,d) + 0.5*weight(j1)*EzRL(i,j,k,j1)*EzyG(j1,d)
                    end do
                    dBx(i,j,k,d) = (dBx(i,j,k,d) - EzVertex(i,j,k)*EzU(d)/hy + EzVertex(i,j - 1,k)*EzD(d)/hy)/mmE(d)
                end do
            end do
        end do
    end do
    
    ! DG scheme of Ez2
    do i = 1,Nx
        do j = 0,Ny
            do k = 0,Nphi
                do d = 1,kk + 1
                    do i1 = 1,NumGLP
                        dBy(i,j,k,d) = dBy(i,j,k,d) - 0.5*weight(i1)*EzUD(i,j,k,i1)*EzxG(i1,d)
                    end do
                    dBy(i,j,k,d) = (dBy(i,j,k,d) + EzVertex(i,j,k)*EzR(d)/hx - EzVertex(i - 1,j,k)*EzL(d)/hx)/mmE(d)
                end do
            end do
        end do
    end do
    
    end subroutine Lh
    
    !*****************************************************************************************************
    
    subroutine HLL_Flux
    
    use com
    
    if (SR < 0) then
        Fhat1 = FR1
    else if (SL > 0) then
        Fhat1 = FL1
    else
        Fhat1 = ( SR*FL1 - SL*FR1 + SL*SR*(UR1 - UL1) )/(SR - SL)
    end if
    
    end subroutine HLL_Flux
    
    !*****************************************************************************************************
    
    subroutine HLLD_Flux
    
    use com
    
    real(8) Sstar,rhoR,rhoL,uRbot,uLbot,PtotR,PtotL,ER,EL,u1R,u1L,u2R,u2L,u3R,u3L,B1R,B1L,B2R,B2L,B3R,B3L
    real(8) rhoRstar,rhoLstar,B1star,B2star,B3star,u1star,u2star,u1Rstar,u1Lstar,u2Rstar,u2Lstar,u3Rstar,u3Lstar
    real(8) BRbot,BLbot,Bstarbot,ptotRstar,ptotLstar,ERstar,ELstar,SRstar,SLstar
    real(8) B1Lstar,B1Rstar,B2Lstar,B2Rstar,B3Lstar,B3Rstar
    real(8) AL,AR,CL,CR,GL,GR,DL,DR,rhoRstarstar,rhoLstarstar
    real(8) u1starstar,u2starstar,u3starstar
    real(8) B1starstar,B2starstar,B3starstar
    real(8) ERstarstar,ELstarstar
    
    if (SR < 0) then
        Fhat1 = FR1
        Ustar = UR1
    else if (SL > 0) then
        Fhat1 = FL1
        Ustar = UL1
    else
        ! Stage 1
        Ustar = ( SR*UR1 - SL*UL1 + FL1 - FR1 )/(SR - SL)
        
        rhoR = UR1(1)
        u1R = UR1(2)/rhoR
        u2R = UR1(3)/rhoR
        u3R = UR1(4)/rhoR
        ER = UR1(5)
        B1R = UR1(6)
        B2R = UR1(7)
        B3R = UR1(8)
        
        rhoL = UL1(1)
        u1L = UL1(2)/rhoL
        u2L = UL1(3)/rhoL
        u3L = UL1(4)/rhoL
        EL = UL1(5)
        B1L = UL1(6)
        B2L = UL1(7)
        B3L = UL1(8)
        
        if (direction == 1) then
            Sstar = Ustar(2)/Ustar(1)
            uRbot = u1R
            uLbot = u1L
            BRbot = B1R
            BLbot = B1L
            Bstarbot = Ustar(6)
        else if (direction == 2) then
            Sstar = Ustar(3)/Ustar(1)
            uRbot = u2R
            uLbot = u2L
            BRbot = B2R
            BLbot = B2L
            Bstarbot = Ustar(7)
        end if
        
        rhoRstar = rhoR*(SR - uRbot)/(SR - Sstar)
        rhoLstar = rhoL*(SL - uLbot)/(SL - Sstar)
        
        SRstar = Sstar + (Bstarbot**2/rhoRstar)**0.5
        SLstar = Sstar - (Bstarbot**2/rhoLstar)**0.5
        
        PtotR = gamma1*(ER - 0.5d0*rhoR*(u1R**2 + u2R**2 + u3R**2))
        PtotL = gamma1*(EL - 0.5d0*rhoL*(u1L**2 + u2L**2 + u3L**2))
        
        PtotRstar = PtotR + rhoR*(SR - uRbot)*(Sstar - uRbot) + Bstarbot**2 - BRbot**2
        PtotLstar = PtotL + rhoL*(SL - uLbot)*(Sstar - uLbot) + Bstarbot**2 - BLbot**2
        
        AR = rhoR*(SR - uRbot)**2 - Bstarbot*BRbot
        AL = rhoL*(SL - uLbot)**2 - Bstarbot*BLbot
        
        CR = rhoR*(SR - uRbot)*(BRbot - Bstarbot)
        CL = rhoL*(SL - uLbot)*(BLbot - Bstarbot)
        
        DR = rhoR*(SR - uRbot)*(SR - Sstar) - Bstarbot**2 + 1e-15
        DL = rhoL*(SL - uLbot)*(SL - Sstar) - Bstarbot**2 + 1e-15
        
        if (direction == 1) then
            B1Rstar = Bstarbot
            B1Lstar = Bstarbot
            B2Rstar = (AR*B2R + CR*u2R)/DR
            B2Lstar = (AL*B2L + CL*u2L)/DL
        else if (direction == 2) then
            B1Rstar = (AR*B1R + CR*u1R)/DR
            B1Lstar = (AL*B1L + CL*u1L)/DL
            B2Rstar = Bstarbot
            B2Lstar = Bstarbot
        end if
        B3Rstar = (AR*B3R + CR*u3R)/DR
        B3Lstar = (AL*B3L + CL*u3L)/DL
        
        if (direction == 1) then
            URstar(2) = (UR1(2)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B1R - Bstarbot*B1Rstar)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot)                     + BRbot*B2R - Bstarbot*B2Rstar)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B1L - Bstarbot*B1Lstar)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot)                     + BLbot*B2L - Bstarbot*B2Lstar)/(SL - Sstar)
        else if (direction == 2) then
            URstar(2) = (UR1(2)*(SR - uRbot)                     + BRbot*B1R - Bstarbot*B1Rstar)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B2R - Bstarbot*B2Rstar)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot)                     + BLbot*B1L - Bstarbot*B1Lstar)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B2L - Bstarbot*B2Lstar)/(SL - Sstar)
        end if
        URstar(4) = (UR1(4)*(SR - uRbot) + BRbot*B3R - Bstarbot*B3Rstar)/(SR - Sstar)
        ULstar(4) = (UL1(4)*(SL - uLbot) + BLbot*B3L - Bstarbot*B3Lstar)/(SL - Sstar)
        
        u1Rstar = URstar(2)/rhoRstar
        u2Rstar = URstar(3)/rhoRstar
        u3Rstar = URstar(4)/rhoRstar
        
        u1Lstar = ULstar(2)/rhoLstar
        u2Lstar = ULstar(3)/rhoLstar
        u3Lstar = ULstar(4)/rhoLstar
        
        ERstar = ( (SR - uRbot)*ER - PtotR*uRbot + PtotRstar*Sstar + BRbot*(u1R*B1R + u2R*B2R + u3R*B3R) - Bstarbot*(u1Rstar*B1Rstar + u2Rstar*B2Rstar + u3Rstar*B3Rstar) )/(SR - Sstar)
        ELstar = ( (SL - uLbot)*EL - PtotL*uLbot + PtotLstar*Sstar + BLbot*(u1L*B1L + u2L*B2L + u3L*B3L) - Bstarbot*(u1Lstar*B1Lstar + u2Lstar*B2Lstar + u3Lstar*B3Lstar) )/(SL - Sstar)
        
        URstar(1) = rhoRstar
        URstar(5) = ERstar
        URstar(6:8) = Ustar(6:8)
        
        ULstar(1) = rhoLstar
        ULstar(5) = ELstar
        ULstar(6:8) = Ustar(6:8)
        
        if ((abs(Bstarbot) > 1e-6) .and. (SRstar - SLstar > 1e-6)) then
            if (SLstar > 0) then
                Fhat1 = FL1 + SL*(ULstar - UL1)
                Ustar = ULstar
            else if (SRstar < 0) then
                Fhat1 = FR1 + SR*(URstar - UR1)
                Ustar = URstar
            else
                ! Stage 2
                FR1 = FR1 + SR*(URstar - UR1)
                FL1 = FL1 + SL*(ULstar - UL1)
            
                rhoRstarstar = rhoRstar
                rhoLstarstar = rhoLstar
            
                u1starstar = (abs(rhoLstar)**0.5*u1Lstar + abs(rhoRstar)**0.5*u1Rstar + sign(1d0,Bstarbot)*(B1Rstar - B1Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                u2starstar = (abs(rhoLstar)**0.5*u2Lstar + abs(rhoRstar)**0.5*u2Rstar + sign(1d0,Bstarbot)*(B2Rstar - B2Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                u3starstar = (abs(rhoLstar)**0.5*u3Lstar + abs(rhoRstar)**0.5*u3Rstar + sign(1d0,Bstarbot)*(B3Rstar - B3Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
            
                B1starstar = (abs(rhoLstar)**0.5*B1Lstar + abs(rhoRstar)**0.5*B1Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u1Rstar - u1Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                B2starstar = (abs(rhoLstar)**0.5*B2Lstar + abs(rhoRstar)**0.5*B2Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u2Rstar - u2Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
                B3starstar = (abs(rhoLstar)**0.5*B3Lstar + abs(rhoRstar)**0.5*B3Rstar + sign(1d0,Bstarbot)*abs(rhoLstar*rhoRstar)**0.5*(u3Rstar - u3Lstar))/(abs(rhoLstar)**0.5 + abs(rhoRstar)**0.5)
            
                ERstarstar = ERstar + abs(rhoRstar)**0.5*(u1Rstar*B1Rstar + u2Rstar*B2Rstar + u3Rstar*B3Rstar - u1starstar*B1starstar - u2starstar*B2starstar - u3starstar*B3starstar)*sign(1d0,Bstarbot)
                ELstarstar = ELstar - abs(rhoLstar)**0.5*(u1Lstar*B1Lstar + u2Lstar*B2Lstar + u3Lstar*B3Lstar - u1starstar*B1starstar - u2starstar*B2starstar - u3starstar*B3starstar)*sign(1d0,Bstarbot)
                
                URstarstar(1) = rhoRstarstar
                URstarstar(2) = rhoRstarstar*u1starstar
                URstarstar(3) = rhoRstarstar*u2starstar
                URstarstar(4) = rhoRstarstar*u3starstar
                URstarstar(5) = ERstarstar
                URstarstar(6) = B1starstar
                URstarstar(7) = B2starstar
                URstarstar(8) = B3starstar
            
                ULstarstar(1) = rhoLstarstar
                ULstarstar(2) = rhoLstarstar*u1starstar
                ULstarstar(3) = rhoLstarstar*u2starstar
                ULstarstar(4) = rhoLstarstar*u3starstar
                ULstarstar(5) = ELstarstar
                ULstarstar(6) = B1starstar
                ULstarstar(7) = B2starstar
                ULstarstar(8) = B3starstar
            
                if (Sstar > 0) then
                    Fhat1 = FL1 + SLstar*(ULstarstar - ULstar)
                    Ustar = ULstarstar
                else
                    Fhat1 = FR1 + SRstar*(URstarstar - URstar)
                    Ustar = URstarstar
                end if
                
            end if
        else
            if (Sstar > 0) then
                Fhat1 = FL1 + SL*(ULstar - UL1)
                Ustar = ULstar
            else
                Fhat1 = FR1 + SR*(URstar - UR1)
                Ustar = URstar
            end if
        end if
        
    end if
    
    end subroutine HLLD_Flux
    
    !*****************************************************************************************************
    
    subroutine HLLD_Flux_2D
    
    use com
    
    real(8) alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real(8) alphayRU,alphayLU,alphayRD,alphayLD
    real(8) alphax2D,alphay2D,SR1,SL1
    real(8) EzRU,EzLU,EzRD,EzLD
    real(8) BxRU,BxLU,BxRD,BxLD
    real(8) ByRU,ByLU,ByRD,ByLD
    real(8) EzR1,EzL1,EzU1,EzD1
    real(8) EzRstar,EzLstar,EzUstar,EzDstar,EzStarStar
    real(8) Ustarstar(NumEq)
    real(8) FRstar(NumEq),FLstar(NumEq),FUstar(NumEq),FDstar(NumEq)
    real(8) FxRU(NumEq),FyRU(NumEq),FxLU(NumEq),FyLU(NumEq),FxRD(NumEq),FyRD(NumEq),FxLD(NumEq),FyLD(NumEq)
    
    call eigenvalueMm(alphaxRU,betaxRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),1,0)
    call eigenvalueMm(alphayRU,betayRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),0,1)
    call eigenvalueMm(alphaxLU,betaxLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),1,0)
    call eigenvalueMm(alphayLU,betayLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),0,1)
    call eigenvalueMm(alphaxRD,betaxRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),1,0)
    call eigenvalueMm(alphayRD,betayRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),0,1)
    call eigenvalueMm(alphaxLD,betaxLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),1,0)
    call eigenvalueMm(alphayLD,betayLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),0,1)
    
    SR = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    SL = min(betaxRU,betaxLU,betaxRD,betaxLD)
    SU = max(alphayRU,alphayLU,alphayRD,alphayLD)
    SD = min(betayRU,betayLU,betayRD,betayLD)
    
    call MHD_flux(URU1,FxRU,FyRU)
    call MHD_flux(ULU1,FxLU,FyLU)
    call MHD_flux(URD1,FxRD,FyRD)
    call MHD_flux(ULD1,FxLD,FyLD)
    
    BxRU = URU1(6)
    BxLU = ULU1(6)
    BxRD = URD1(6)
    BxLD = ULD1(6)
    
    ByRU = URU1(7)
    ByLU = ULU1(7)
    ByRD = URD1(7)
    ByLD = ULD1(7)
    
    direction = 1
    
    UR1 = URU1
    UL1 = ULU1
    FR1 = FxRU
    FL1 = FxLU
    call HLLD_Flux
    FUstar = Fhat1
    UUstar = Ustar
    
    UR1 = URD1
    UL1 = ULD1
    FR1 = FxRD
    FL1 = FxLD
    call HLLD_Flux
    FDstar = Fhat1
    UDstar = Ustar
    
    SR1 = SR
    SL1 = SL
    SR = SU
    SL = SD
    direction = 2
    
    UR1 = URU1
    UL1 = URD1
    FR1 = FyRU
    FL1 = FyRD
    call HLLD_Flux
    FRstar = Fhat1
    URstar = Ustar
    
    UR1 = ULU1
    UL1 = ULD1
    FR1 = FyLU
    FL1 = FyLD
    call HLLD_Flux
    FLstar = Fhat1
    ULstar = Ustar
    
    SR = SR1
    SL = SL1
    
    BxUstar = UUstar(6)
    BxDstar = UDstar(6)
    ByRstar = URstar(7)
    ByLstar = ULstar(7)
    
    EzRstar = FRstar(6)
    EzLstar = FLstar(6)
    EzUstar = -FUstar(7)
    EzDstar = -FDstar(7)
    
    EzStarStar = 0.25*(EzRstar + EzLstar + EzUstar + EzDstar)
    
    if (SL > 0) then
        Ezhat = EzLstar
    else if (SR < 0) then
        Ezhat = EzRstar
    else if (SD > 0) then
        Ezhat = EzDstar
    else if (SU < 0) then
        Ezhat = EzUstar
    else
        Ezhat = EzStarStar
    end if
    
    end subroutine HLLD_Flux_2D
    
    !*****************************************************************************************************
    
    subroutine HLL_Flux_2D
    
    use com
    
    real(8) alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real(8) alphayRU,alphayLU,alphayRD,alphayLD
    real(8) alphax2D,alphay2D
    real(8) EzRU,EzLU,EzRD,EzLD
    real(8) BxRU,BxLU,BxRD,BxLD
    real(8) ByRU,ByLU,ByRD,ByLD
    real(8) EzR1,EzL1,EzU1,EzD1
    real(8) EzRstar,EzLstar,EzUstar,EzDstar,EzStarStar
    real(8) Ustarstar(NumEq)
    real(8) FRstar(NumEq),FLstar(NumEq),FUstar(NumEq),FDstar(NumEq)
    real(8) FxRU(NumEq),FyRU(NumEq),FxLU(NumEq),FyLU(NumEq),FxRD(NumEq),FyRD(NumEq),FxLD(NumEq),FyLD(NumEq)
    
    call eigenvalueMm(alphaxRU,betaxRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),1,0)
    call eigenvalueMm(alphayRU,betayRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),0,1)
    call eigenvalueMm(alphaxLU,betaxLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),1,0)
    call eigenvalueMm(alphayLU,betayLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),0,1)
    call eigenvalueMm(alphaxRD,betaxRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),1,0)
    call eigenvalueMm(alphayRD,betayRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),0,1)
    call eigenvalueMm(alphaxLD,betaxLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),1,0)
    call eigenvalueMm(alphayLD,betayLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),0,1)
    
    SR = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    SL = min(betaxRU,betaxLU,betaxRD,betaxLD)
    SU = max(alphayRU,alphayLU,alphayRD,alphayLD)
    SD = min(betayRU,betayLU,betayRD,betayLD)
    
    call MHD_flux(URU1,FxRU,FyRU)
    call MHD_flux(ULU1,FxLU,FyLU)
    call MHD_flux(URD1,FxRD,FyRD)
    call MHD_flux(ULD1,FxLD,FyLD)
    
    BxRU = URU1(6)
    BxLU = ULU1(6)
    BxRD = URD1(6)
    BxLD = ULD1(6)
    
    ByRU = URU1(7)
    ByLU = ULU1(7)
    ByRD = URD1(7)
    ByLD = ULD1(7)
    
    call HLL1(URU1,ULU1,FxRU,FxLU,SR,SL,UUstar,FUstar)
    call HLL1(URD1,ULD1,FxRD,FxLD,SR,SL,UDstar,FDstar)
    call HLL1(URU1,URD1,FyRU,FyRD,SU,SD,URstar,FRstar)
    call HLL1(ULU1,ULD1,FyLU,FyLD,SU,SD,ULstar,FLstar)
    
    BxUstar = UUstar(6)
    BxDstar = UDstar(6)
    ByRstar = URstar(7)
    ByLstar = ULstar(7)
    
    EzRstar = FRstar(6)
    EzLstar = FLstar(6)
    EzUstar = -FUstar(7)
    EzDstar = -FDstar(7)
    
    BxStarStar = ( 2*SR*SU*BxRU - 2*SL*SU*BxLU + 2*SL*SD*BxLD - 2*SR*SD*BxRD - SR*(EzRU - EzRD) + SL*(EzLU - EzLD) - (SR - SL)*(EzUstar - EzDstar) )/( 2*(SR - SL)*(SU - SD) )
    ByStarStar = ( 2*SR*SU*ByRU - 2*SL*SU*ByLU + 2*SL*SD*ByLD - 2*SR*SD*ByRD + SU*(EzRU - EzLU) - SD*(EzRD - EzLD) + (SU - SD)*(EzRstar - EzLstar) )/( 2*(SR - SL)*(SU - SD) )
    
    EzStarStar = 0.25*(EzRstar + EzLstar + EzUstar + EzDstar) - 0.25*SU*(BxUstar - BxStarStar) - 0.25*SD*(BxDstar - BxStarStar) + 0.25*SR*(ByRstar - ByStarStar) + 0.25*SL*(ByLstar - ByStarStar)
    
    if (SL > 0) then
        Ezhat = EzLstar
    else if (SR < 0) then
        Ezhat = EzRstar
    else if (SD > 0) then
        Ezhat = EzDstar
    else if (SU < 0) then
        Ezhat = EzUstar
    else
        Ezhat = EzStarStar
    end if
    
    end subroutine HLL_Flux_2D
    
    !*****************************************************************************************************
    
    subroutine HLL1(UR,UL,FR,FL,SR,SL,Uhat,Fhat)
    
    real(8) UR(8),UL(8),FR(8),FL(8),SR,SL,Uhat(8),Fhat(8)
    
    if (SR < 0) then
        Fhat = FR
        Uhat = UR
    else if (SL > 0) then
        Fhat = FL
        Uhat = UL
    else
        Fhat = ( SR*FL - SL*FR + SL*SR*(UR - UL) )/(SR - SL)
        Uhat = ( SR*UR - SL*UL - (FR - FL) )/(SR - SL)
    end if
    
    end subroutine HLL1
    
    !*****************************************************************************************************
    
    subroutine MHD_flux(Uh,Fx,Fy)
    
    real(8) Uh(8),Fx(8),Fy(8)
    real(8) rho,u,v,w,p,E,B1,B2,B3,gamma,gamma1
    real(8) S,T,K
      
    gamma = 5d0/3d0
    gamma1 = gamma - 1
    
    rho = Uh(1)
    u = Uh(2)/rho
    v = Uh(3)/rho
    w = Uh(4)/rho
    E = Uh(5)
    B1 = Uh(6)
    B2 = Uh(7)
    B3 = Uh(8)

    p = gamma1*(E - 0.5*rho*(u**2 + v**2 + w**2) - 0.5*(B1**2 + B2**2 + B3**2))
    S = p + 0.5*(B1**2 + B2**2 + B3**2)
    T = E + S
    K = u*B1 + v*B2 + w*B3
    
    Fx(1) = Uh(2)
    Fx(2) = rho*u**2 + S - B1**2
    Fx(3) = rho*u*v - B1*B2
    Fx(4) = rho*u*w - B1*B3
    Fx(5) = T*u - K*B1
    Fx(6) = 0
    Fx(7) = u*B2 - v*B1
    Fx(8) = u*B3 - w*B1
    
    Fy(1) = Uh(3)
    Fy(2) = rho*u*v - B1*B2
    Fy(3) = rho*v**2 + S - B2**2
    Fy(4) = rho*w*v - B3*B2
    Fy(5) = T*v - K*B2
    Fy(6) = v*B1 - u*B2
    Fy(7) = 0
    Fy(8) = v*B3 - w*B2
    
    end subroutine MHD_flux
    
    !*****************************************************************************************************
    
    subroutine HLLC_Flux
    
    use com
    
    real(8) Sstar,rhoR,rhoL,uRbot,uLbot,PtotR,PtotL,ER,EL,u1R,u1L,u2R,u2L,u3R,u3L,B1R,B1L,B2R,B2L,B3R,B3L
    real(8) rhoRstar,rhoLstar,B1star,B2star,B3star,u1star,u2star,u1Rstar,u1Lstar,u2Rstar,u2Lstar,u3Rstar,u3Lstar
    real(8) BRbot,BLbot,Bstarbot,Ptotstar,ERstar,ELstar
    
    if (SR < 0) then
        Fhat1 = FR1
        Ustar = UR1
    else if (SL > 0) then
        Fhat1 = FL1
        Ustar = UL1
    else
        Ustar = ( SR*UR1 - SL*UL1 + FL1 - FR1 )/(SR - SL)
        
        rhoR = UR1(1)
        u1R = UR1(2)/rhoR
        u2R = UR1(3)/rhoR
        u3R = UR1(4)/rhoR
        ER = UR1(5)
        B1R = UR1(6)
        B2R = UR1(7)
        B3R = UR1(8)
        
        rhoL = UL1(1)
        u1L = UL1(2)/rhoL
        u2L = UL1(3)/rhoL
        u3L = UL1(4)/rhoL
        EL = UL1(5)
        B1L = UL1(6)
        B2L = UL1(7)
        B3L = UL1(8)
        
        B1star = Ustar(6)
        B2star = Ustar(7)
        B3star = Ustar(8)
        
        if (direction == 1) then
            Sstar = Ustar(2)/Ustar(1)
            uRbot = u1R
            uLbot = u1L
            BRbot = B1R
            BLbot = B1L
            Bstarbot = B1star
        else if (direction == 2) then
            Sstar = Ustar(3)/Ustar(1)
            uRbot = u2R
            uLbot = u2L
            BRbot = B2R
            BLbot = B2L
            Bstarbot = B2star
        end if
        
        rhoRstar = rhoR*(SR - uRbot)/(SR - Sstar)
        rhoLstar = rhoL*(SL - uLbot)/(SL - Sstar)
        
        PtotR = gamma1*(ER - 0.5d0*rhoR*(u1R**2 + u2R**2 + u3R**2))
        PtotL = gamma1*(EL - 0.5d0*rhoL*(u1L**2 + u2L**2 + u3L**2))
        
        PtotRstar = PtotR + rhoR*(SR - uRbot)*(Sstar - uRbot) + Bstarbot**2 - BRbot**2
        PtotLstar = PtotL + rhoL*(SL - uLbot)*(Sstar - uLbot) + Bstarbot**2 - BLbot**2
        
        if (direction == 1) then
            URstar(2) = (UR1(2)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B1R - Bstarbot*B1star)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot)                     + BRbot*B2R - Bstarbot*B2star)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B1L - Bstarbot*B1star)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot)                     + BLbot*B2L - Bstarbot*B2star)/(SL - Sstar)
        else if (direction == 2) then
            URstar(2) = (UR1(2)*(SR - uRbot)                     + BRbot*B1R - Bstarbot*B1star)/(SR - Sstar)
            URstar(3) = (UR1(3)*(SR - uRbot) + ptotRstar - ptotR + BRbot*B2R - Bstarbot*B2star)/(SR - Sstar)
            ULstar(2) = (UL1(2)*(SL - uLbot)                     + BLbot*B1L - Bstarbot*B1star)/(SL - Sstar)
            ULstar(3) = (UL1(3)*(SL - uLbot) + ptotLstar - ptotL + BLbot*B2L - Bstarbot*B2star)/(SL - Sstar)
        end if
        URstar(4) = (UR1(4)*(SR - uRbot) + BRbot*B3R - Bstarbot*B3star)/(SR - Sstar)
        ULstar(4) = (UL1(4)*(SL - uLbot) + BLbot*B3L - Bstarbot*B3star)/(SL - Sstar)
        
        u1Rstar = URstar(2)/rhoRstar
        u2Rstar = URstar(3)/rhoRstar
        u3Rstar = URstar(4)/rhoRstar
        
        u1Lstar = ULstar(2)/rhoLstar
        u2Lstar = ULstar(3)/rhoLstar
        u3Lstar = ULstar(4)/rhoLstar
        
        ERstar = ( (SR - uRbot)*ER - PtotR*uRbot + PtotRstar*Sstar + BRbot*(u1R*B1R + u2R*B2R + u3R*B3R) - Bstarbot*(u1Rstar*B1star + u2Rstar*B2star + u3Rstar*B3star) )/(SR - Sstar)
        ELstar = ( (SL - uLbot)*EL - PtotL*uLbot + PtotLstar*Sstar + BLbot*(u1L*B1L + u2L*B2L + u3L*B3L) - Bstarbot*(u1Lstar*B1star + u2Lstar*B2star + u3Lstar*B3star) )/(SL - Sstar)
        
        URstar(1) = rhoRstar
        URstar(5) = ERstar
        URstar(6:8) = Ustar(6:8)
        
        ULstar(1) = rhoLstar
        ULstar(5) = ELstar
        ULstar(6:8) = Ustar(6:8)
    
        if (Sstar > 0) then
            Fhat1 = FL1 + SL*(ULstar - UL1)
            Ustar = ULstar
        else
            Fhat1 = FR1 + SR*(URstar - UR1)
            Ustar = URstar
        end if
        
    end if
    
    end subroutine HLLC_Flux
    
    !*****************************************************************************************************
    
    subroutine HLLC_Flux_2D
    
    use com
    
    real(8) alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real(8) alphayRU,alphayLU,alphayRD,alphayLD
    real(8) alphax2D,alphay2D,SR1,SL1
    real(8) EzRU,EzLU,EzRD,EzLD
    real(8) BxRU,BxLU,BxRD,BxLD
    real(8) ByRU,ByLU,ByRD,ByLD
    real(8) EzR1,EzL1,EzU1,EzD1
    real(8) EzRstar,EzLstar,EzUstar,EzDstar,EzStarStar
    real(8) Ustarstar(NumEq)
    real(8) FRstar(NumEq),FLstar(NumEq),FUstar(NumEq),FDstar(NumEq)
    real(8) FxRU(NumEq),FyRU(NumEq),FxLU(NumEq),FyLU(NumEq),FxRD(NumEq),FyRD(NumEq),FxLD(NumEq),FyLD(NumEq)
    
    call eigenvalueMm(alphaxRU,betaxRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),1,0)
    call eigenvalueMm(alphayRU,betayRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),0,1)
    call eigenvalueMm(alphaxLU,betaxLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),1,0)
    call eigenvalueMm(alphayLU,betayLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),0,1)
    call eigenvalueMm(alphaxRD,betaxRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),1,0)
    call eigenvalueMm(alphayRD,betayRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),0,1)
    call eigenvalueMm(alphaxLD,betaxLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),1,0)
    call eigenvalueMm(alphayLD,betayLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),0,1)
    
    SR = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    SL = min(betaxRU,betaxLU,betaxRD,betaxLD)
    SU = max(alphayRU,alphayLU,alphayRD,alphayLD)
    SD = min(betayRU,betayLU,betayRD,betayLD)
    
    call MHD_flux(URU1,FxRU,FyRU)
    call MHD_flux(ULU1,FxLU,FyLU)
    call MHD_flux(URD1,FxRD,FyRD)
    call MHD_flux(ULD1,FxLD,FyLD)
    
    BxRU = URU1(6)
    BxLU = ULU1(6)
    BxRD = URD1(6)
    BxLD = ULD1(6)
    
    ByRU = URU1(7)
    ByLU = ULU1(7)
    ByRD = URD1(7)
    ByLD = ULD1(7)
    
    direction = 1
    
    UR1 = URU1
    UL1 = ULU1
    FR1 = FxRU
    FL1 = FxLU
    call HLLC_Flux
    FUstar = Fhat1
    UUstar = Ustar
    
    UR1 = URD1
    UL1 = ULD1
    FR1 = FxRD
    FL1 = FxLD
    call HLLC_Flux
    FDstar = Fhat1
    UDstar = Ustar
    
    SR1 = SR
    SL1 = SL
    SR = SU
    SL = SD
    direction = 2
    
    UR1 = URU1
    UL1 = URD1
    FR1 = FyRU
    FL1 = FyRD
    call HLLC_Flux
    FRstar = Fhat1
    URstar = Ustar
    
    UR1 = ULU1
    UL1 = ULD1
    FR1 = FyLU
    FL1 = FyLD
    call HLLC_Flux
    FLstar = Fhat1
    ULstar = Ustar
    
    SR = SR1
    SL = SL1
    
    BxUstar = UUstar(6)
    BxDstar = UDstar(6)
    ByRstar = URstar(7)
    ByLstar = ULstar(7)
    
    EzRstar = FRstar(6)
    EzLstar = FLstar(6)
    EzUstar = -FUstar(7)
    EzDstar = -FDstar(7)
    
    EzStarStar = 0.25*(EzRstar + EzLstar + EzUstar + EzDstar)
    
    if (SL > 0) then
        Ezhat = EzLstar
    else if (SR < 0) then
        Ezhat = EzRstar
    else if (SD > 0) then
        Ezhat = EzDstar
    else if (SU < 0) then
        Ezhat = EzUstar
    else
        Ezhat = EzStarStar
    end if
    
    end subroutine HLLC_Flux_2D
    
    !*****************************************************************************************************
    
    subroutine LF_Flux
    
    use com
    
    Fhat1 = 0.5d0*(FR1 + FL1 - max(abs(SR),abs(SL))*(UR1 - UL1))
    
    end subroutine LF_Flux
    
    !*****************************************************************************************************
    
    subroutine LF_Flux_2D
    
    use com
    
    real(8) alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real(8) alphayRU,alphayLU,alphayRD,alphayLD
    real(8) alphax2D,alphay2D
    real(8) EzRU,EzLU,EzRD,EzLD
    real(8) B1RU,B1LU,B1RD,B1LD
    real(8) B2RU,B2LU,B2RD,B2LD
    real(8) EzR1,EzL1,EzU1,EzD1
    
    alphax2D = max(abs(URU1(2)/URU1(1)),abs(ULU1(2)/ULU1(1)),abs(URD1(2)/URD1(1)),abs(ULD1(2)/ULD1(1)))
    alphay2D = max(abs(URU1(3)/URU1(1)),abs(ULU1(3)/ULU1(1)),abs(URD1(3)/URD1(1)),abs(ULD1(3)/ULD1(1)))
    
    call calculate_Ez(EzRU,URU1(2)/URU1(1),URU1(3)/URU1(1),URU1(6),URU1(7))
    call calculate_Ez(EzLU,ULU1(2)/ULU1(1),ULU1(3)/ULU1(1),ULU1(6),ULU1(7))
    call calculate_Ez(EzRD,URD1(2)/URD1(1),URD1(3)/URD1(1),URD1(6),URD1(7))
    call calculate_Ez(EzLD,ULD1(2)/ULD1(1),ULD1(3)/ULD1(1),ULD1(6),ULD1(7))
    
    B1RU = URU1(6)
    B1LU = ULU1(6)
    B1RD = URD1(6)
    B1LD = ULD1(6)
    
    B2RU = URU1(7)
    B2LU = ULU1(7)
    B2RD = URD1(7)
    B2LD = ULD1(7)
    
    Ezhat = 0.25*(EzRU + EzLU + EzRD + EzLD) - 0.25*alphay2D*(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)) + 0.25*alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    
    end subroutine LF_Flux_2D
    
    !*****************************************************************************************************
    
    subroutine calculate_Ez(Ez,u1,u2,B1,B2)
    
    real(8) Ez,u1,u2,B1,B2
    
    Ez = u2*B1 - u1*B2
    
    end subroutine calculate_Ez
    
    !*****************************************************************************************************
    
    subroutine div_free_Balsara
    
    use com
    
    real(8) BxR(kk + 1),BxL(kk + 1),ByU(kk + 1),ByD(kk + 1)
    real(8) a0R,a1R,a2R,a0L,a1L,a2L
    real(8) b0U,b1U,b2U,b0D,b1D,b2D
    real(8) a00,a10,a01,a20,a11,a02,a30,a12
    real(8) b00,b10,b01,b20,b11,b02,b21,b03
    real(8) rxy,ryx
    
    uh(:,:,:,:,6:7) = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
            
                BxR = Bx(i,j,k,:)
                BxL = Bx(i - 1,j,k,:)
                ByU = By(i,j,k,:)
                ByD = By(i,j - 1,k,:)
            
                a0R = BxR(1)
                a1R = BxR(2)
                a2R = BxR(3)
    
                a0L = BxL(1)
                a1L = BxL(2)
                a2L = BxL(3)
    
                b0U = ByU(1)
                b1U = ByU(2)
                b2U = ByU(3)
    
                b0D = ByD(1)
                b1D = ByD(2)
                b2D = ByD(3)
            
                rxy = hx/hy
                ryx = hy/hx
    
                ! The reconstruction of B = (Bx,By) from the interface
                a00 = (1d0/2d0)*(a0R + a0L) + (1d0/6d0)*rxy*(b1U - b1D)
                a10 = (1d0/2d0)*(a0R - a0L) + (1d0/15d0)*rxy*(b2U - b2D)
                a01 = (1d0/2d0)*(a1R + a1L)
                a20 = -(1d0/4d0)*rxy*(b1U - b1D)
                a11 = (1d0/2d0)*(a1R - a1L)
                a02 = (1d0/2d0)*(a2R + a2L)
                a30 = -(1d0/6d0)*rxy*(b2U - b2D)
                a12 = (1d0/2d0)*(a2R - a2L)
            
                b00 = (1d0/2d0)*(b0U + b0D) + (1d0/6d0)*ryx*(a1R - a1L)
                b01 = (1d0/2d0)*(b0U - b0D) + (1d0/15d0)*ryx*(a2R - a2L)
                b10 = (1d0/2d0)*(b1U + b1D)
                b02 = -(1d0/4d0)*ryx*(a1R - a1L)
                b11 = (1d0/2d0)*(b1U - b1D)
                b20 = (1d0/2d0)*(b2U + b2D)
                b03 = -(1d0/6d0)*ryx*(a2R - a2L)
                b21 = (1d0/2d0)*(b2U - b2D)
            
                uh(i,j,k,1,6) = a00
                uh(i,j,k,2,6) = a10
                uh(i,j,k,3,6) = a01
                uh(i,j,k,4,6) = a20
                uh(i,j,k,5,6) = a11
                uh(i,j,k,6,6) = a02
                uh(i,j,k,7,6) = a30
                uh(i,j,k,9,6) = a12
            
                uh(i,j,k,1,7) = b00
                uh(i,j,k,2,7) = b10
                uh(i,j,k,3,7) = b01
                uh(i,j,k,4,7) = b20
                uh(i,j,k,5,7) = b11
                uh(i,j,k,6,7) = b02
                uh(i,j,k,8,7) = b21
                uh(i,j,k,10,7) = b03
            
            end do
        end do
    end do
    
    end subroutine div_free_Balsara
    
    !*****************************************************************************************************