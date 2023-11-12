    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP,dimPk
    parameter(N_process = 16)
    parameter(Nx0 = 100, Ny0 = 100, Lphi = 0, kk = 3, NumEq = 8, NumGLP = 5, RKorder = 4, flux_type = 2)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) pi,ly,gamma,gamma1
    parameter(dimPk = (kk + 2)*(kk + 3)/2)
    parameter(dimPk1 = (kk + 1)*(kk + 2)/2)
    parameter(Nphi = max(2*Lphi - 1,0))
    parameter(Nphi1 = Nphi + 1)
    !parameter(gamma = 1.4d0) ! other
    parameter(gamma = 5d0/3d0) ! jet, R-T
    parameter(gamma1 = gamma - 1)
    parameter(pi = 4*atan(1d0))
    parameter(ly = 32d0/pi)
    
    ! The numerical solution and mesh
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,tRK,t1,t2,alphax,alphay,totaldiv,rij,eta,nu
    real(8) uh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),du(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uhs(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) phs(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),qhs(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) p1s(Nx,Ny,0:Nphi,dimPk,NumEq),p2s(Nx,Ny,0:Nphi,dimPk,NumEq)
    real(8) q1s(Nx,Ny,0:Nphi,dimPk,NumEq),q2s(Nx,Ny,0:Nphi,dimPk,NumEq)
    real(8) uinit(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uI(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),uII(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uh00(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) hx,hy,Xc(Nx),Yc(Ny),Xc0(Nx0),Yc0(Ny0),Phi(0:Nphi),hx1,hy1,hphi
    real(8) Bx(0:Nx,0:Ny1,0:Nphi,kk + 1),By(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) dBx(0:Nx,0:Ny1,0:Nphi,kk + 1),dBy(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxI(0:Nx,0:Ny1,0:Nphi,kk + 1),ByI(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxII(0:Nx,0:Ny1,0:Nphi,kk + 1),ByII(0:Nx1,0:Ny,0:Nphi,kk + 1)
    
    ! The basis
    real(8) lambda(NumGLP),weight(NumGLP),sink(0:Nphi,Lphi),cosk(0:Nphi,Lphi)
    real(8) phiG(NumGLP,NumGLP,dimPk),phixG(NumGLP,NumGLP,dimPk),phiyG(NumGLP,NumGLP,dimPk),mm(dimPk)
    real(8) phiGLL(NumGLP,NumGLP,dimPk,2),lambdaL(NumGLP)
    real(8) phiGR(NumGLP,dimPk), phiGL(NumGLP,dimPk), phiGU(NumGLP,dimPk), phiGD(NumGLP,dimPk)
    real(8) phiRU(dimPk), phiLU(dimPk), phiRD(dimPk), phiLD(dimPk)
    real(8) EzG(NumGLP,kk + 1),EzxG(NumGLP,kk + 1),EzyG(NumGLP,kk + 1),mmE(kk + 1)
    real(8) EzR(kk + 1),EzL(kk + 1),EzU(kk + 1),EzD(kk + 1),omega1(Nx,Ny,0:Nphi)
    real(8) EzRL(0:Nx,Ny,0:Nphi,NumGLP), EzUD(Nx,0:Ny,0:Nphi,NumGLP), EzVertex(0:Nx,0:Ny,0:Nphi)
    
    ! The Lh
    real(8) uGint3D(NumGLP,NumGLP,0:Nphi,NumEq),uGint(NumGLP,NumGLP,NumEq)
    real(8) pGint3D(NumGLP,NumGLP,0:Nphi,NumEq),qGint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) mu1(NumGLP,NumGLP,0:Nphi,NumEq),mu2(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHSC(NumGLP,NumGLP,0:Nphi,NumEq),RHSCopen,RG(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHS(NumGLP,NumGLP,0:Nphi,NumEq),Fzsin(Lphi),Fzcos(Lphi),Fzzsin(Lphi),Fzzcos(Lphi)
    real(8) FR1(NumEq),FL1(NumEq),UR1(NumEq),UL1(NumEq),Fhat1(NumEq),SR,SL,SU,SD
    real(8) URstar(NumEq),ULstar(NumEq),Ustar(NumEq),UUstar(NumEq),UDstar(NumEq)
    real(8) URU1(NumEq),ULU1(NumEq),URD1(NumEq),ULD1(NumEq)
    real(8) URUs1(NumEq),ULUs1(NumEq),URDs1(NumEq),ULDs1(NumEq)
    real(8) PRU1(NumEq),PLU1(NumEq),PRD1(NumEq),PLD1(NumEq)
    real(8) QRU1(NumEq),QLU1(NumEq),QRD1(NumEq),QLD1(NumEq)
    real(8) URstarstar(NumEq),ULstarstar(NumEq),Ezhat,Hhat1(NumEq)
    real(8) URU(0:Nx1,0:Ny1,0:Nphi,NumEq),ULU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URD(0:Nx1,0:Ny1,0:Nphi,NumEq),ULD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URUs(0:Nx1,0:Ny1,0:Nphi,NumEq),ULUs(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URDs(0:Nx1,0:Ny1,0:Nphi,NumEq),ULDs(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) PRU(0:Nx1,0:Ny1,0:Nphi,NumEq),PLU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) PRD(0:Nx1,0:Ny1,0:Nphi,NumEq),PLD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) QRU(0:Nx1,0:Ny1,0:Nphi,NumEq),QLU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) QRD(0:Nx1,0:Ny1,0:Nphi,NumEq),QLD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) L2(NumEq),L2pre(NumEq)
    real(8) Hij(NumEq),p1ij(NumEq),p2ij(NumEq),q1ij(NumEq),q2ij(NumEq)
    real(8) p1sij(NumEq),p2sij(NumEq),q1sij(NumEq),q2sij(NumEq),uhsij(NumEq)
    real(8) p1wij(NumEq),p2wij(NumEq),q1wij(NumEq),q2wij(NumEq)
    real(8) HRU(NumEq),HLU(NumEq),HRD(NumEq),HLD(NumEq)
    real(8) p1sGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq),p2sGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq)
    real(8) q1sGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq),q2sGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq)
    real(8) psGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq),qsGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq)
    real(8) usGint3D(NumGLP,NumGLP,Nx,Ny,0:Nphi,NumEq)
    real(8) URs(0:Nx,Ny,0:Nphi,NumGLP,NumEq)
    real(8) ULs(Nx1,Ny,0:Nphi,NumGLP,NumEq)
    real(8) UUs(Nx,0:Ny,0:Nphi,NumGLP,NumEq)
    real(8) UDs(Nx,Ny1,0:Nphi,NumGLP,NumEq)
    
    ! The Limiter
    real(8) M,beta
    real(8) DeltaUR1(NumEq,1),DeltaUL1(NumEq,1),DeltaUU1(NumEq,1),DeltaUD1(NumEq,1),DeltaU1(NumEq,1),DeltaUmod1(NumEq,1)
    real(8) DeltaUR1mod(NumEq,1),DeltaUL1mod(NumEq,1),DeltaUU1mod(NumEq,1),DeltaUD1mod(NumEq,1)
    real(8) R(NumEq,NumEq),L(NumEq,NumEq)
    real(8) DeltaUR(NumEq,1),DeltaUL(NumEq,1),DeltaU(NumEq,1),DeltaUmod(NumEq,1)
    real(8) Is_trouble_cell(Nx,Ny,0:Nphi)
    
    integer bcR,bcL,bcU,bcD,direction
    integer myid,myid1,the_id,the_id2
    integer myidx,myidy,the_idx,the_idy
    integer numprocs, namelen, rc,ierr,status(MPI_STATUS_SIZE),myid0
    character * (MPI_MAX_PROCESSOR_NAME) processor_name
    
    end module com
    
    !*****************************************************************************************************
    
    ! Smooth vortex
    module init1
    
    use com
    
    contains
    
    function r2(x,y)
    real(8) r2,x,y
    r2 = x**2 + y**2
    end function r2
    
    function p(x,y,z)
    real(8) x,y,p,z
    p = 1 - r2(x,y)/(8.0d0*pi**2)*exp(1 - r2(x,y))
    end function p
    
    function rho(x,y,z)
    real(8) x,y,rho,z
    rho = 1
    end function rho
    
    function v1(x,y,z)
    real(8) v1,x,y,z
    v1 = 1 - 1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y)))
    end function v1
    
    function v2(x,y,z)
    real(8) v2,x,y,z
    v2 = 1 + 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y)))
    end function v2
    
    function v3(x,y,z)
    real(8) v3,x,y,z
    v3 = 0
    end function v3
    
    function B1(x,y,z)
    real(8) B1,x,y,z
    B1 = -1d0/(2.0d0*pi)*y*exp(0.5d0*(1 - r2(x,y)))
    end function B1
    
    function B2(x,y,z)
    real(8) B2,x,y,z
    B2 = 1d0/(2.0d0*pi)*x*exp(0.5d0*(1 - r2(x,y)))
    end function B2
    
    function B3(x,y,z)
    real(8) B3,x,y,z
    B3 = 0
    end function B3
    
    ! Us
    
    function ps(x,y,z)
    real(8) x,y,ps,z
    ps = 0
    end function ps
    
    function rhos(x,y,z)
    real(8) x,y,rhos,z
    rhos = 0
    end function rhos
    
    function v1s(x,y,z)
    real(8) v1s,x,y,z
    v1s = 0
    end function v1s
    
    function v2s(x,y,z)
    real(8) v2s,x,y,z
    v2s = 0
    end function v2s
    
    function v3s(x,y,z)
    real(8) v3s,x,y,z
    v3s = 0
    end function v3s
    
    function B1s(x,y,z)
    real(8) B1s,x,y,z
    B1s = 0
    end function B1s
    
    function B2s(x,y,z)
    real(8) B2s,x,y,z
    B2s = 0
    end function B2s
    
    function B3s(x,y,z)
    real(8) B3s,x,y,z
    B3s = 0
    end function B3s
    
    
    
    subroutine mesh
    
    use com
    
    eta = 0
    nu = 0.01
    
    xa = -10
    xb = 10
    ya = -10
    yb = 10

    bcL = 1
    bcR = 1
    bcU = 1
    bcD = 1
    
    RHSCopen = 1

    tend = 2
    
    end subroutine mesh
    
    end module init1
    
    !*****************************************************************************************************
    
    ! sin
    module init2
    
    use com
    
    contains
    
    function p(x,y,z)
    real(8) x,y,p,z
    p = 1
    end function p
    
    function rho(x,y,z)
    real(8) x,y,rho,z
    rho = 1 + 0.5*sin(x + y)
    end function rho
    
    function v1(x,y,z)
    real(8) v1,x,y,z
    v1 = 1
    end function v1
    
    function v2(x,y,z)
    real(8) v2,x,y,z
    v2 = 1
    end function v2
    
    function v3(x,y,z)
    real(8) v3,x,y,z
    v3 = 0
    end function v3
    
    function B1(x,y,z)
    real(8) B1,x,y,z
    B1 = -0.5*sin(y)
    end function B1
    
    function B2(x,y,z)
    real(8) B2,x,y,z
    B2 = 0.5*sin(x)
    end function B2
    
    function B3(x,y,z)
    real(8) B3,x,y,z
    B3 = 0
    end function B3
    
    ! Us
    
    function ps(x,y,z)
    real(8) x,y,ps,z
    ps = 0
    end function ps
    
    function rhos(x,y,z)
    real(8) x,y,rhos,z
    rhos = 0
    end function rhos
    
    function v1s(x,y,z)
    real(8) v1s,x,y,z
    v1s = 0
    end function v1s
    
    function v2s(x,y,z)
    real(8) v2s,x,y,z
    v2s = 0
    end function v2s
    
    function v3s(x,y,z)
    real(8) v3s,x,y,z
    v3s = 0
    end function v3s
    
    function B1s(x,y,z)
    real(8) B1s,x,y,z
    B1s = 0
    end function B1s
    
    function B2s(x,y,z)
    real(8) B2s,x,y,z
    B2s = 0
    end function B2s
    
    function B3s(x,y,z)
    real(8) B3s,x,y,z
    B3s = 0
    end function B3s
    
    
    
    subroutine mesh
    
    use com
    
    eta = 0.05
    nu = 0.01
    
    xa = 0
    xb = 2*pi
    ya = 0
    yb = 2*pi

    bcL = 1
    bcR = 1
    bcU = 1
    bcD = 1
    
    RHSCopen = 1

    tend = 5
    
    end subroutine mesh
    
    end module init2
    
    !*****************************************************************************************************
    
    ! Orszag-Tang vortex
    module init3
    
    use com
    
    contains
    
    function p(x,y,z)
    real(8) x,y,p,z
    p = gamma
    end function p
    
    function rho(x,y,z)
    real(8) x,y,rho,z
    rho = gamma**2
    end function rho
    
    function v1(x,y,z)
    real(8) v1,x,y,z
    v1 = -sin(y)
    end function v1
    
    function v2(x,y,z)
    real(8) v2,x,y,z
    v2 = sin(x)
    end function v2
    
    function v3(x,y,z)
    real(8) v3,x,y,z
    v3 = 0
    end function v3
    
    function B1(x,y,z)
    real(8) B1,x,y,z
    B1 = -sin(y)
    end function B1
    
    function B2(x,y,z)
    real(8) B2,x,y,z
    B2 = sin(2*x)
    end function B2
    
    function B3(x,y,z)
    real(8) B3,x,y,z
    B3 = 0
    end function B3
    
    subroutine mesh
    
    use com
    
    eta = 0.05
    nu = 0.01
    
    xa = 0
    xb = 2*pi
    ya = 0
    yb = 2*pi

    bcL = 1
    bcR = 1
    bcU = 1
    bcD = 1

    tend = 0.00067
    
    end subroutine mesh
    
    end module init3
    
    !*****************************************************************************************************
    
    ! 2D-Tearing
    module init4
    
    use com
    
    contains
    
    function p(x,y,z)
    real(8) x,y,p,z
    p = 0
    end function p
    
    function rho(x,y,z)
    real(8) x,y,rho,z
    rho = 0
    end function rho
    
    function v1(x,y,z)
    real(8) v1,x,y,z
    v1 = 0
    end function v1
    
    function v2(x,y,z)
    real(8) v2,x,y,z
    v2 = 0
    end function v2
    
    function v3(x,y,z)
    real(8) v3,x,y,z
    v3 = 0
    end function v3
    
    function B1(x,y,z)
    real(8) B1,x,y,z
    B1 = -0.01/ly*exp(-(x/3d0)**2)*cos(y/ly)
    end function B1
    
    function B2(x,y,z)
    real(8) B2,x,y,z
    B2 = - 0.01*2*x/(3d0**2)*exp(-(x/3d0)**2)*sin(y/ly)
    end function B2
    
    function B3(x,y,z)
    real(8) B3,x,y,z
    B3 = 0
    end function B3
    
    ! Us
    
    function v1s(x,y,z)
    real(8) v1s,x,y,z
    v1s = 0
    end function v1s
    
    function v2s(x,y,z)
    real(8) v2s,x,y,z
    v2s = 0
    end function v2s
    
    function v3s(x,y,z)
    real(8) v3s,x,y,z
    v3s = 0
    end function v3s
    
    function B1s(x,y,z)
    real(8) B1s,x,y,z
    B1s = 0
    end function B1s
    
    function B2s(x,y,z)
    real(8) x,y,z,B2s
    B2s = tanh(x/3d0)
    end function B2s
    
    function B3s(x,y,z)
    real(8) B3s,x,y,z
    B3s = 0
    end function B3s
    
    function ps(x,y,z)
    real(8) x,y,ps,z
    ps = 0.25 + 0.5*(1 - B2s(x,y,z)**2)
    end function ps
    
    function rhos(x,y,z)
    real(8) x,y,rhos,z
    rhos = 4*ps(x,y,z)
    end function rhos
    
    subroutine mesh
    
    use com
    
    eta = 0.05
    nu = 0.01
    
    xa = -32
    xb = 32
    ya = 0
    yb = 64

    bcL = 2
    bcR = 2
    bcU = 1
    bcD = 1

    tend = 1000
    
    end subroutine mesh
    
    end module init4
    
    !*****************************************************************************************************
    
    ! 2D-Tearing
    module init5
    
    use com
    
    contains
    
    function p(x,y,z)
    real(8) x,y,p,z
    p = 0
    end function p
    
    function rho(x,y,z)
    real(8) x,y,rho,z
    rho = 0
    end function rho
    
    function v1(x,y,z)
    real(8) v1,x,y,z
    v1 = 0
    end function v1
    
    function v2(x,y,z)
    real(8) v2,x,y,z
    v2 = 0
    end function v2
    
    function v3(x,y,z)
    real(8) v3,x,y,z
    v3 = 0
    end function v3
    
    function B1(x,y,z)
    real(8) B1,x,y,z
    B1 = 0
    end function B1
    
    function B2(x,y,z)
    real(8) B2,x,y,z
    B2 = 0
    end function B2
    
    function B3(x,y,z)
    real(8) B3,x,y,z
    B3 = 0
    end function B3
    
    ! Us
    
    function v1s(x,y,z)
    real(8) v1s,x,y,z
    v1s = 0
    end function v1s
    
    function v2s(x,y,z)
    real(8) v2s,x,y,z
    v2s = 0
    end function v2s
    
    function v3s(x,y,z)
    real(8) v3s,x,y,z
    v3s = 0
    end function v3s
    
    function B1s(x,y,z)
    real(8) B1s,x,y,z
    B1s = 0
    end function B1s
    
    function B2s(x,y,z)
    real(8) x,y,z,B2s
    B2s = tanh(x/3d0)
    end function B2s
    
    function B3s(x,y,z)
    real(8) B3s,x,y,z
    B3s = 0
    end function B3s
    
    function ps(x,y,z)
    real(8) x,y,ps,z
    ps = 0.25 + 0.5*(1 - B2s(x,y,z)**2)
    end function ps
    
    function rhos(x,y,z)
    real(8) x,y,rhos,z
    rhos = 4*ps(x,y,z)
    end function rhos
    
    subroutine mesh
    
    use com
    
    eta = 0
    nu = 0
    
    xa = -32
    xb = 32
    ya = 0
    yb = 64

    bcL = 2
    bcR = 2
    bcU = 1
    bcD = 1

    tend = 10
    
    end subroutine mesh
    
    end module init5
    
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
    
    print *,"process",myid1,"is alive,the index is",myidx,myidy
    
    call init_data
    
    if (RKorder == 1) then
        call Euler_Forward
    else if (RKorder == 3) then
        call RK3
    else if (RKorder == 4) then
        call RK4
    end if
    
    call set_bc
    
    call calculate_L2_Error
    
    call save_solution
    
    call save_solution2
    
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
    
    ! The initial value:
    ! 1: smooth vortex
    ! 2: sin
    ! 3: Orszag-Tang vortex
    ! 4: 2D-Tearing
    ! 5: 2D-Tearing s
    
    use init4
    
    real(8) U1
    U1(x,y,z) = rho(x,y,z) + rhos(x,y,z)
    real(8) U2
    U2(x,y,z) = v1(x,y,z) + v1s(x,y,z)
    real(8) U3
    U3(x,y,z) = v2(x,y,z) + v2s(x,y,z)
    real(8) U4
    U4(x,y,z) = v3(x,y,z) + v3s(x,y,z)
    real(8) U5
    U5(x,y,z) = p(x,y,z) + ps(x,y,z)
    real(8) U6
    U6(x,y,z) = B1(x,y,z) + B1s(x,y,z)
    real(8) U7
    U7(x,y,z) = B2(x,y,z) + B2s(x,y,z)
    real(8) U8
    U8(x,y,z) = B3(x,y,z) + B3s(x,y,z)
    
    real(8) U1s
    U1s(x,y,z) = rhos(x,y,z)
    real(8) U2s
    U2s(x,y,z) = v1s(x,y,z)
    real(8) U3s
    U3s(x,y,z) = v2s(x,y,z)
    real(8) U4s
    U4s(x,y,z) = v3s(x,y,z)
    real(8) U5s
    U5s(x,y,z) = ps(x,y,z)
    real(8) U6s
    U6s(x,y,z) = B1s(x,y,z)
    real(8) U7s
    U7s(x,y,z) = B2s(x,y,z)
    real(8) U8s
    U8s(x,y,z) = B3s(x,y,z)
    
    call mesh
    
    hx = (xb - xa)/Nx0
    hy = (yb - ya)/Ny0
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    hphi = 2d0*pi/Nphi1
    
    do i = 1,Nx0
        Xc0(i) = xa + (i - 0.5)*hx
    end do
    Xc = Xc0((myidx - 1)*Nx + 1:myidx*Nx)
    
    do j = 1,Ny0
        Yc0(j) = ya + (j - 0.5)*hy
    end do
    Yc = Yc0((myidy - 1)*Ny + 1:myidy*Ny)
    
    do k = 0,Nphi
        Phi(k) = k*hphi
    end do
    
    call get_basis
    
    uh = 0
    uinit = 0
    uhs = 0
    
    ! L2 Pro for Uh
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh(i,j,k,d,1) = uh(i,j,k,d,1) + 0.25*weight(i1)*weight(j1)*U1(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,2) = uh(i,j,k,d,2) + 0.25*weight(i1)*weight(j1)*U2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,3) = uh(i,j,k,d,3) + 0.25*weight(i1)*weight(j1)*U3(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,4) = uh(i,j,k,d,4) + 0.25*weight(i1)*weight(j1)*U4(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,5) = uh(i,j,k,d,5) + 0.25*weight(i1)*weight(j1)*U5(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,6) = uh(i,j,k,d,6) + 0.25*weight(i1)*weight(j1)*U6(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,7) = uh(i,j,k,d,7) + 0.25*weight(i1)*weight(j1)*U7(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,8) = uh(i,j,k,d,8) + 0.25*weight(i1)*weight(j1)*U8(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            
                            uhs(i,j,k,d,1) = uhs(i,j,k,d,1) + 0.25*weight(i1)*weight(j1)*U1s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,2) = uhs(i,j,k,d,2) + 0.25*weight(i1)*weight(j1)*U2s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,3) = uhs(i,j,k,d,3) + 0.25*weight(i1)*weight(j1)*U3s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,4) = uhs(i,j,k,d,4) + 0.25*weight(i1)*weight(j1)*U4s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,5) = uhs(i,j,k,d,5) + 0.25*weight(i1)*weight(j1)*U5s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,6) = uhs(i,j,k,d,6) + 0.25*weight(i1)*weight(j1)*U6s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,7) = uhs(i,j,k,d,7) + 0.25*weight(i1)*weight(j1)*U7s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uhs(i,j,k,d,8) = uhs(i,j,k,d,8) + 0.25*weight(i1)*weight(j1)*U8s(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        uh(:,:,:,d,:) = uh(:,:,:,d,:)/mm(d)
        uhs(:,:,:,d,:) = uhs(:,:,:,d,:)/mm(d)
    end do
    
    ! L2 Pro for Bx and By
    Bx0 = 0
    By0 = 0
    do i = 0,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,kk + 1
                    do j1 = 1,NumGLP
                        Bx(i,j,k,d) = Bx(i,j,k,d) + 0.5*weight(j1)*EzG(j1,d)*U6(xa + ((myidx - 1)*Nx + i)*hx,Yc(j) + hy1*lambda(j1),Phi(k))
                    end do
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do k = 0,Nphi
                do d = 1,kk + 1
                    do i1 = 1,NumGLP
                        By(i,j,k,d) = By(i,j,k,d) + 0.5*weight(i1)*EzG(i1,d)*U7(Xc(i) + hx1*lambda(i1),ya + ((myidy - 1)*Ny + j)*hy,Phi(k))
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,kk + 1
        Bx(:,:,:,d) = Bx(:,:,:,d)/mmE(d)
        By(:,:,:,d) = By(:,:,:,d)/mmE(d)
    end do
    
    call calculate_pqs
    
    end subroutine init_data
    
    !*****************************************************************************************************
    
    subroutine save_solution
    
    use com
    
    real uhsave(NumEq)
    integer the_idx1,the_idy1
    
    !uh = du
    !uh = phs
    
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
        
        do i = 1,Nx0
            write(9,*) Xc0(i)
        end do
        
        do j = 1,Ny0
            write(10,*) Yc0(j)
        end do
        
        close(9)
        close(10)
        
    end if
    
    do d = 1,1
        do j = 1,Ny0
            do i = 1,Nx0
                    
                the_idx1 = mod(i,Nx)
                if (the_idx1 == 0) then
                    the_idx1 = Nx
                end if
                the_idx = (i - the_idx1)/Nx + 1
            
                the_idy1 = mod(j,Ny)
                if (the_idy1 == 0) then
                    the_idy1 = Ny
                end if
                the_idy = (j - the_idy1)/Ny + 1
            
                the_id = the_idx + Nx_process*(the_idy - 1)
                    
                if (the_id /= 1) then
                    if (myid1 == the_id) then
                        call MPI_SEND(uh(the_idx1,the_idy1,0,d,:),NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
                    end if
                    if (myid1 == 1) then
                        call MPI_RECV(uhsave,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
                    end if
                else if (the_id == 1) then
                    if (myid1 == 1) then
                        uhsave = uh(the_idx1,the_idy1,0,d,:)
                    end if
                end if
            
                if (myid1 == 1) then
                    write(1,*) uhsave(1)
                    write(2,*) uhsave(2)
                    write(3,*) uhsave(3)
                    write(4,*) uhsave(4)
                    write(5,*) uhsave(5)
                    write(6,*) uhsave(6)
                    write(7,*) uhsave(7)
                    write(8,*) uhsave(8)
                end if
                    
            end do
        end do
    end do
    
    if (myid1 == 1) then
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
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
        lambda(4) = 0.8611363115940525752239465d0
        
        weight(1) = 0.3478548451374538573730639d0
        weight(2) = 0.6521451548625461426269361d0
        weight(3) = 0.6521451548625461426269361d0   
        weight(4) = 0.3478548451374538573730639d0
    else if (NumGLP == 5) then
        lambda(1) = -0.9061798459386639927976269d0     
        lambda(2) = -0.5384693101056830910363144d0     
        lambda(3) = 0d0                                 
        lambda(4) = 0.5384693101056830910363144d0     
        lambda(5) = 0.9061798459386639927976269d0     
        
        weight(1) = 0.2369268850561890875142640d0
        weight(2) = 0.4786286704993664680412915d0
        weight(3) = 0.5688888888888888888888889d0
        weight(4) = 0.4786286704993664680412915d0
        weight(5) = 0.2369268850561890875142640d0
        
        lambdaL(1) = -1
        lambdaL(2) = -0.6546536707079771437983
        lambdaL(3) = 0
        lambdaL(4) = 0.654653670707977143798
        lambdaL(5) = 1
    else if (NumGLP == 6) then
        lambda(1) = -0.9324695142031520278123016d0     
        lambda(2) = -0.6612093864662645136613996d0    
        lambda(3) = -0.2386191860831969086305017d0     
        !lambda(4) = 0.2386191860831969086305017d0     
        !lambda(5) = 0.6612093864662645136613996d0     
        !lambda(6) = 0.9324695142031520278123016d0     
        
        weight(1) = 0.1713244923791703450402961d0
        weight(2) = 0.3607615730481386075698335d0
        weight(3) = 0.4679139345726910473898703d0
        !weight(4) = 0.4679139345726910473898703d0
        !weight(5) = 0.3607615730481386075698335d0
        !weight(6) = 0.1713244923791703450402961d0
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
    phiRU(11) = 8d0/35d0
    phiRU(12) = 2d0/5d0
    phiRU(13) = 4d0/9d0
    phiRU(14) = 2d0/5d0
    phiRU(15) = 8d0/35d0
    
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
    phiLU(11) = 8d0/35d0
    phiLU(12) = -2d0/5d0
    phiLU(13) = 4d0/9d0
    phiLU(14) = -2d0/5d0
    phiLU(15) = 8d0/35d0
    
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
    phiRD(11) = 8d0/35d0
    phiRD(12) = -2d0/5d0
    phiRD(13) = 4d0/9d0
    phiRD(14) = -2d0/5d0
    phiRD(15) = 8d0/35d0
    
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
    phiLD(11) = 8d0/35d0
    phiLD(12) = 2d0/5d0
    phiLD(13) = 4d0/9d0
    phiLD(14) = 2d0/5d0
    phiLD(15) = 8d0/35d0
    
    do i = 1,NumGLP
        do j = 1,NumGLP
            phiG(i,j,1) = 1
            phiGLL(i,j,1,1) = 1
            phiGLL(i,j,1,2) = 1
            phiGR(j,1) = 1
            phiGL(j,1) = 1
            phiGU(i,1) = 1
            phiGD(i,1) = 1
            phixG(i,j,1) = 0
            phiyG(i,j,1) = 0
            mm(1) = 1
            
            phiG(i,j,2) = lambda(i)
            phiGLL(i,j,2,1) = lambdaL(i)
            phiGLL(i,j,2,2) = lambda(i)
            phiGR(j,2) = 1
            phiGL(j,2) = -1
            phiGU(i,2) = lambda(i)
            phiGD(i,2) = lambda(i)
            phixG(i,j,2) = 1d0/hx1
            phiyG(i,j,2) = 0
            mm(2) = 1d0/3d0
                
            phiG(i,j,3) = lambda(j)
            phiGLL(i,j,3,1) = lambda(j)
            phiGLL(i,j,3,2) = lambdaL(j)
            phiGR(j,3) = lambda(j)
            phiGL(j,3) = lambda(j)
            phiGU(i,3) = 1
            phiGD(i,3) = -1
            phixG(i,j,3) = 0
            phiyG(i,j,3) = 1d0/hy1
            mm(3) = 1d0/3d0
                
            phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
            phiGLL(i,j,4,1) = lambdaL(i)**2 - 1d0/3d0
            phiGLL(i,j,4,2) = lambda(i)**2 - 1d0/3d0
            phiGR(j,4) = 2d0/3d0
            phiGL(j,4) = 2d0/3d0
            phiGU(i,4) = lambda(i)**2 - 1d0/3d0
            phiGD(i,4) = lambda(i)**2 - 1d0/3d0
            phixG(i,j,4) = 2d0*lambda(i)/hx1
            phiyG(i,j,4) = 0
            mm(4) = 4d0/45d0
                    
            phiG(i,j,5) = lambda(i)*lambda(j)
            phiGLL(i,j,5,1) = lambdaL(i)*lambda(j)
            phiGLL(i,j,5,2) = lambda(i)*lambdaL(j)
            phiGR(j,5) = lambda(j)
            phiGL(j,5) = -lambda(j)
            phiGU(i,5) = lambda(i)
            phiGD(i,5) = -lambda(i)
            phixG(i,j,5) = lambda(j)/hx1
            phiyG(i,j,5) = lambda(i)/hy1
            mm(5) = 1d0/9d0
                    
            phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,1) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,2) = lambdaL(j)**2 - 1d0/3d0
            phiGR(j,6) = lambda(j)**2 - 1d0/3d0
            phiGL(j,6) = lambda(j)**2 - 1d0/3d0
            phiGU(i,6) = 2d0/3d0
            phiGD(i,6) = 2d0/3d0
            phixG(i,j,6) = 0
            phiyG(i,j,6) = 2d0*lambda(j)/hy1
            mm(6) = 4d0/45d0
                
            phiG(i,j,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGLL(i,j,7,1) = lambdaL(i)**3 - 3d0*lambdaL(i)/5d0
            phiGLL(i,j,7,2) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGR(j,7) = 2d0/5d0
            phiGL(j,7) = -2d0/5d0
            phiGU(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGD(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phixG(i,j,7) = (3*lambda(i)**2 - 3d0/5d0)/hx1
            phiyG(i,j,7) = 0
            mm(7) = 4d0/175d0
                
            phiG(i,j,8) = (lambda(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,1) = (lambdaL(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,2) = (lambda(i)**2 - 1d0/3d0)*(lambdaL(j))
            phiGR(j,8) = (2d0/3d0)*(lambda(j))
            phiGL(j,8) = (2d0/3d0)*(lambda(j))
            phiGU(i,8) = (lambda(i)**2 - 1d0/3d0)
            phiGD(i,8) = -(lambda(i)**2 - 1d0/3d0)
            phixG(i,j,8) = 2d0*lambda(i)*lambda(j)/hx1
            phiyG(i,j,8) = (lambda(i)**2 - 1d0/3d0)/hy1
            mm(8) = 4d0/135d0
                
            phiG(i,j,9) = (lambda(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,1) = (lambdaL(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,2) = (lambda(i))*(lambdaL(j)**2 - 1d0/3d0)
            phiGR(j,9) = (lambda(j)**2 - 1d0/3d0)
            phiGL(j,9) = -(lambda(j)**2 - 1d0/3d0)
            phiGU(i,9) = lambda(i)*(2d0/3d0)
            phiGD(i,9) = lambda(i)*(2d0/3d0)
            phixG(i,j,9) = (lambda(j)**2 - 1d0/3d0)/hx1
            phiyG(i,j,9) = 2d0*lambda(i)*lambda(j)/hy1
            mm(9) = 4d0/135d0
                
            phiG(i,j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,1) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,2) = lambdaL(j)**3 - 3d0*lambdaL(j)/5d0
            phiGR(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGL(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGU(i,10) = 2d0/5d0
            phiGD(i,10) = -2d0/5d0
            phixG(i,j,10) = 0
            phiyG(i,j,10) = (3*lambda(j)**2 - 3d0/5d0)/hy1
            mm(10) = 4d0/175d0
            
            phiG(i,j,11) = lambda(i)**4 - 6*lambda(i)**2/7d0 + 3d0/35d0
            phiGLL(i,j,11,1) = lambdaL(i)**4 - 6*lambdaL(i)**2/7d0 + 3d0/35d0
            phiGLL(i,j,11,2) = lambda(i)**4 - 6*lambda(i)**2/7d0 + 3d0/35d0
            phiGR(j,11) = 8d0/35d0
            phiGL(j,11) = 8d0/35d0
            phiGU(i,11) = lambda(i)**4 - 6*lambda(i)**2/7d0 + 3d0/35d0
            phiGD(i,11) = lambda(i)**4 - 6*lambda(i)**2/7d0 + 3d0/35d0
            phixG(i,j,11) = (4*lambda(i)**3 - 12*lambda(i)/7d0)/hx1
            phiyG(i,j,11) = 0
            mm(11) = 64d0/11025d0
            
            phiG(i,j,12) = (lambda(i)**3 - 3d0*lambda(i)/5d0)*lambda(j)
            phiGLL(i,j,12,1) = (lambdaL(i)**3 - 3d0*lambdaL(i)/5d0)*lambda(j)
            phiGLL(i,j,12,2) = (lambda(i)**3 - 3d0*lambda(i)/5d0)*lambdaL(j)
            phiGR(j,12) = (2d0/5d0)*lambda(j)
            phiGL(j,12) = -(2d0/5d0)*lambda(j)
            phiGU(i,12) = (lambda(i)**3 - 3d0*lambda(i)/5d0)
            phiGD(i,12) = -(lambda(i)**3 - 3d0*lambda(i)/5d0)
            phixG(i,j,12) = (3*lambda(i)**2 - 3d0/5d0)*lambda(j)/hx1
            phiyG(i,j,12) = (lambda(i)**3 - 3d0*lambda(i)/5d0)/hy1
            mm(12) = 4d0/525d0
            
            phiG(i,j,13) = (lambda(i)**2 - 1d0/3d0)*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,13,1) = (lambdaL(i)**2 - 1d0/3d0)*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,13,2) = (lambda(i)**2 - 1d0/3d0)*(lambdaL(j)**2 - 1d0/3d0)
            phiGR(j,13) = (2d0/3d0)*(lambda(j)**2 - 1d0/3d0)
            phiGL(j,13) = (2d0/3d0)*(lambda(j)**2 - 1d0/3d0)
            phiGU(i,13) = (lambda(i)**2 - 1d0/3d0)*(2d0/3d0)
            phiGD(i,13) = (lambda(i)**2 - 1d0/3d0)*(2d0/3d0)
            phixG(i,j,13) = 2*lambda(i)*(lambda(j)**2 - 1d0/3d0)/hx1
            phiyG(i,j,13) = (lambda(i)**2 - 1d0/3d0)*2*lambda(j)/hy1
            mm(13) = 16d0/2025d0
            
            phiG(i,j,14) = lambda(i)*(lambda(j)**3 - 3d0*lambda(j)/5d0)
            phiGLL(i,j,14,1) = lambdaL(i)*(lambda(j)**3 - 3d0*lambda(j)/5d0)
            phiGLL(i,j,14,2) = lambda(i)*(lambdaL(j)**3 - 3d0*lambdaL(j)/5d0)
            phiGR(j,14) = (lambda(j)**3 - 3d0*lambda(j)/5d0)
            phiGL(j,14) = -(lambda(j)**3 - 3d0*lambda(j)/5d0)
            phiGU(i,14) = lambda(i)*(2d0/5d0)
            phiGD(i,14) = -lambda(i)*(2d0/5d0)
            phixG(i,j,14) = (lambda(j)**3 - 3d0*lambda(j)/5d0)/hx1
            phiyG(i,j,14) = lambda(i)*(3*lambda(j)**2 - 3d0/5d0)/hy1
            mm(14) = 4d0/525d0
            
            phiG(i,j,15) = lambda(j)**4 - 6*lambda(j)**2/7d0 + 3d0/35d0
            phiGLL(i,j,15,1) = lambda(j)**4 - 6*lambda(j)**2/7d0 + 3d0/35d0
            phiGLL(i,j,15,2) = lambdaL(j)**4 - 6*lambdaL(j)**2/7d0 + 3d0/35d0
            phiGR(j,15) = lambda(j)**4 - 6*lambda(j)**2/7d0 + 3d0/35d0
            phiGL(j,15) = lambda(j)**4 - 6*lambda(j)**2/7d0 + 3d0/35d0
            phiGU(j,15) = 8d0/35d0
            phiGD(j,15) = 8d0/35d0
            phixG(i,j,15) = 0
            phiyG(i,j,15) = (4*lambda(j)**3 - 12*lambda(j)/7d0)/hy1
            mm(15) = 64d0/11025d0
        end do
        
        EzG(i,1) = 1
        EzG(i,2) = lambda(i)
        EzG(i,3) = lambda(i)**2 - 1d0/3d0
        EzG(i,4) = lambda(i)**3 - 3d0/5d0*lambda(i)
        
        EzxG(i,1) = 0
        EzxG(i,2) = 1/hx1
        EzxG(i,3) = 2*lambda(i)/hx1
        EzxG(i,4) = (3*lambda(i)**2 - 3d0/5d0)/hx1
        
        EzyG(i,1) = 0
        EzyG(i,2) = 1/hy1
        EzyG(i,3) = 2*lambda(i)/hy1
        EzyG(i,4) = (3*lambda(i)**2 - 3d0/5d0)/hy1
        
    end do
    
    EzR(1) = 1
    EzR(2) = 1
    EzR(3) = 2d0/3d0
    EzR(4) = 2d0/5d0
    
    EzL(1) = 1
    EzL(2) = -1
    EzL(3) = 2d0/3d0
    EzL(4) = -2d0/5d0
    
    EzU = EzR
    EzD = EzL
    
    mmE(1) = 1
    mmE(2) = 1d0/3d0
    mmE(3) = 4d0/45d0
    mmE(4) = 4d0/175d0
    
    end subroutine get_basis
    
    !*****************************************************************************************************
    
    subroutine calculate_L2_Error
    
    use com
    
    use init1
    
    real(8) U1
    U1(x,y,z) = rho(x,y,z)
    real(8) U2
    U2(x,y,z) = v1(x,y,z)
    real(8) U3
    U3(x,y,z) = v2(x,y,z)
    real(8) U4
    U4(x,y,z) = v3(x,y,z)
    real(8) U5
    U5(x,y,z) = p(x,y,z)
    real(8) U6
    U6(x,y,z) = B1(x,y,z)
    real(8) U7
    U7(x,y,z) = B2(x,y,z)
    real(8) U8
    U8(x,y,z) = B3(x,y,z)
    
    L2 = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                
                uGint = 0
                do d = 1,dimPk
                    do n = 1,NumEq
                        uGint(:,:,n) = uGint(:,:,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
                
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        L2(1) = L2(1) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,1) - U1(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(2) = L2(2) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,2) - U2(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(3) = L2(3) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,3) - U3(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(4) = L2(4) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,4) - U4(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(5) = L2(5) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,5) - U5(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(6) = L2(6) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,6) - exp(-eta*tend)*U6(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(7) = L2(7) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,7) - exp(-eta*tend)*U7(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(8) = L2(8) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,8) - U8(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                    end do
                end do
                
            end do
        end do
    end do
    
    do the_id = 2,N_process
        
        if (myid1 == the_id) then
            call MPI_SEND(L2,NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(L2pre,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            L2 = L2 + L2pre
        end if
        
    end do
    
    if (myid1 == 1) then
        L2 = (L2/(Nx0*Ny0*Nphi1))**0.5d0
    end if
    
    if (myid1 == 1) then
        print *,"The L2 Error:"
        print *,"rho    : ",L2(1)
        print *,"ux     : ",L2(2)
        print *,"uy     : ",L2(3)
        print *,"uz     : ",L2(4)
        print *,"p      : ",L2(5)
        print *,"Bx     : ",L2(6)
        print *,"By     : ",L2(7)
        print *,"Bz     : ",L2(8)
    end if
    
    end subroutine calculate_L2_Error
    
    !*****************************************************************************************************
    
    subroutine Euler_Forward
    
    use com
    
    CFL = 0.01
    t = 0
    
    call calculate_umax
    
    if (myid1 == 1) then
        print *,t,umax
    end if
    
    do while (t < tend)
        
        call calculate_dt
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        tRK = t
        call Lh
        
        uh = uh + dt*du
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
        end if
        
    end do
    
    end subroutine Euler_Forward
    
    !*****************************************************************************************************
    
    subroutine RK3
    
    use com
    
    CFL = 0.1
    t = 0
    
    do while (t < tend)
        
        call calculate_dt
        
        tRK = t
        
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
        tRK = tRK + dt
        call Lh
        
        uII = (3d0/4d0)*uh00 + (1d0/4d0)*uh + (1d0/4d0)*dt*du
        
        uh = uII
        
        ! Stage III
        tRK = tRK - 0.5*dt
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
    
    CFL = 0.6
    t = 0
    
    call div_free_Balsara
    
    call calculate_umax
    
    call calculate_totaldiv
        
    if (myid1 == 1) then
        open(unit = 12,file = 'Latest_result.txt')
        open(unit = 13,file = 'divT.txt')
        open(unit = 14,file = 'maxB.txt')
        print *,t,umax,totaldiv
        write(12,*) t,umax,totaldiv
        write(13,*) totaldiv
        write(14,*) umax
    end if
    
    do while (t < tend)
        
        call calculate_dt
        tRK = t
        
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
            tRK = tRK + (dt/6d0)
            
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
        
        tRK = tRK - 0.5*dt
        
        ! Stage II
        do i = 6,9
            
            call Lh
            
            uI = uh + (dt/6d0)*du
            BxI = Bx + (dt/6d0)*dBx
            ByI = By + (dt/6d0)*dBy
            
            tRK = tRK + dt/6d0
            
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
        
        call calculate_totaldiv
        
        if (myid1 == 1) then
            print *,t,umax,totaldiv
            write(12,*) t,umax,totaldiv
            write(13,*) totaldiv
            write(14,*) umax
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
                if (abs(uh(i,j,k,1,6)) > umax) then
                    umax = abs(uh(i,j,k,1,6))
                end if
            end do
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
    
    dt = CFL/(alphax/hx + alphay/hy + eta/hx**2 + eta/hy**2)
    
    end subroutine calculate_dt
    
    !*****************************************************************************************************
    
    subroutine set_bc
    
    use com
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny_process
                    do i = 1,Nx_process
                        
                        the_id = i + Nx_process*(j - 1)
                        
                        ! The Uh
                        ! The Right condition
                        if (i == Nx_process) then
                            
                            if (bcR == 1) then ! periodic
                                the_idx = 1
                                the_idy = j
                                the_id2 = the_idx + Nx_process*(the_idy - 1)
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                                end if
                                if (myid1 == the_id) then
                                    call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                                end if
                            else if (bcR == 2) then ! reflection
                                if (myid1 == the_id) then
                                    uh(Nx1,:,k,d,n) = uh(Nx,:,k,d,n)
                                    uh(Nx1,:,k,d,2) = 0
                                    uh(Nx1,:,k,d,6) = 0
                                    !do jj = 0,Ny1
                                    !    call evenex_x(uh(Nx1,jj,k,:,1),uh(Nx,jj,k,:,1))
                                    !    call evenex_x(uh(Nx1,jj,k,:,2),uh(Nx,jj,k,:,2))
                                    !    call evenex_x(uh(Nx1,jj,k,:,3),uh(Nx,jj,k,:,3))
                                    !    call evenex_x(uh(Nx1,jj,k,:,4),uh(Nx,jj,k,:,4))
                                    !    call evenex_x(uh(Nx1,jj,k,:,5),uh(Nx,jj,k,:,5))
                                    !    call evenex_x(uh(Nx1,jj,k,:,6),uh(Nx,jj,k,:,6))
                                    !    call evenex_x(uh(Nx1,jj,k,:,7),uh(Nx,jj,k,:,7))
                                    !    call evenex_x(uh(Nx1,jj,k,:,8),uh(Nx,jj,k,:,8))
                                    !end do
                                end if
                            end if
            
                        else
                
                            the_idx = i + 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if
            
                        
                        
                        ! The Left condition
                        if (i == 1) then
                            
                            if (bcL == 1) then
                                the_idx = Nx_process
                                the_idy = j
                                the_id2 = the_idx + Nx_process*(the_idy - 1)
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                                end if
                                if (myid1 == the_id) then
                                    call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                                end if
                            else if (bcL == 2) then ! reflection
                                if (myid1 == the_id) then
                                    uh(0,:,k,d,n) = uh(1,:,k,d,n)
                                    uh(0,:,k,d,2) = 0
                                    uh(0,:,k,d,6) = 0
                                    !do jj = 0,Ny1
                                    !    call evenex_x(uh(0,jj,k,:,1),uh(1,jj,k,:,1))
                                    !    call evenex_x(uh(0,jj,k,:,2),uh(1,jj,k,:,2))
                                    !    call evenex_x(uh(0,jj,k,:,3),uh(1,jj,k,:,3))
                                    !    call evenex_x(uh(0,jj,k,:,4),uh(1,jj,k,:,4))
                                    !    call evenex_x(uh(0,jj,k,:,5),uh(1,jj,k,:,5))
                                    !    call evenex_x(uh(0,jj,k,:,6),uh(1,jj,k,:,6))
                                    !    call evenex_x(uh(0,jj,k,:,7),uh(1,jj,k,:,7))
                                    !    call evenex_x(uh(0,jj,k,:,8),uh(1,jj,k,:,8))
                                    !end do
                                end if
                            else if (bcL == 5) then
                                if (myid1 == the_id) then
                                    do jj = 0,Ny1
                                        call evenex_x(uh(0,jj,k,:,1),uh(1,jj,k,:,1))
                                        call oddex_x(uh(0,jj,k,:,2),uh(1,jj,k,:,2))
                                        call evenex_x(uh(0,jj,k,:,3),uh(1,jj,k,:,3))
                                        call evenex_x(uh(0,jj,k,:,4),uh(1,jj,k,:,4))
                                        call evenex_x(uh(0,jj,k,:,5),uh(1,jj,k,:,5))
                                        call evenex_x(uh(0,jj,k,:,6),uh(1,jj,k,:,6))
                                        call oddex_x(uh(0,jj,k,:,7),uh(1,jj,k,:,7))
                                        call evenex_x(uh(0,jj,k,:,8),uh(1,jj,k,:,8))
                                    end do
                                end if
                            end if
                        else
                            the_idx = i - 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if
            
                        
            
                        ! The Up condition
                        if (j == Ny_process) then
                            
                            if (bcU == 1) then
                                the_idx = i
                                the_idy = 1
                                the_id2 = the_idx + Nx_process*(the_idy - 1)
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                                end if
                                if (myid1 == the_id) then
                                    call MPI_RECV(uh(0:Nx1,Ny1,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                                end if
                            else if (bcD == 2) then ! reflection
                                if (myid1 == the_id) then
                                    uh(:,Ny1,k,d,n) = uh(:,Ny,k,d,n)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j + 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(0:Nx1,Ny1,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                            end if
                
                        end if
            
                        
            
                        ! The Down condition
                        if (j == 1) then
                
                            if (bcD == 1) then
                                the_idx = i
                                the_idy = Ny_process
                                the_id2 = the_idx + Nx_process*(the_idy - 1)
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                                end if
                                if (myid1 == the_id) then
                                    call MPI_RECV(uh(0:Nx1,0,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                                end if
                            else if (bcD == 2) then ! reflection
                                if (myid1 == the_id) then
                                    uh(:,0,k,d,n) = uh(:,1,k,d,n)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j - 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(0:Nx1,0,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                            end if
                
                        end if
                        
                        
                        
                    end do
                end do
            end do
        end do 
    end do
    
    do d = 1,kk + 1
        do k = 0,Nphi
            do j = 1,Ny_process
                do i = 1,Nx_process
                    
                    the_id = i + Nx_process*(j - 1)
                        
                    ! The Bx
                    ! The Up condition
                    if (j == Ny_process) then
                        
                        if (bcU == 1) then
                            the_idx = i
                            the_idy = 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(Bx(0:Nx,1,k,d),Nx + 1,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(Bx(0:Nx,Ny1,k,d),Nx + 1,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                            end if
                        else if (bcU == 2) then
                            if (myid1 == the_id) then
                                Bx(0:Nx,Ny1,k,d) = Bx(0:Nx,Ny,k,d)
                            end if
                        end if
                    else
                        the_idx = i
                        the_idy = j + 1
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
        
                        if (myid1 == the_id2) then
                            call MPI_SEND(Bx(0:Nx,1,k,d),Nx + 1,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(Bx(0:Nx,Ny1,k,d),Nx + 1,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                        end if
                
                    end if
                    
                    ! The Down condition
                    if (j == 1) then
                
                        if (bcD == 1) then
                            the_idx = i
                            the_idy = Ny_process
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(Bx(0:Nx,Ny,k,d),Nx + 1,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(Bx(0:Nx,0,k,d),Nx + 1,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                            end if
                        else if (bcD == 2) then
                            if (myid1 == the_id) then
                                Bx(0:Nx,0,k,d) = Bx(0:Nx,1,k,d)
                            end if
                        end if
                    else
                        the_idx = i
                        the_idy = j - 1
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                        if (myid1 == the_id2) then
                            call MPI_SEND(Bx(0:Nx,Ny,k,d),Nx + 1,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(Bx(0:Nx,0,k,d),Nx + 1,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                        end if
                
                    end if
                        
                    ! The By
                    ! The Right condition
                    if (i == Nx_process) then
                            
                
                        if (bcR == 1) then
                            the_idx = 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(By(1,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(By(Nx1,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                            end if
                        else if (bcR == 2) then
                            if (myid1 == the_id) then
                                By(Nx1,0:Ny,k,d) = By(Nx,0:Ny,k,d)
                            end if
                        else if (bcR == 4) then
                            the_idx = 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(By(1,0:Ny,k,d)*xb/xa,Ny + 1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(By(Nx1,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if
            
                    else
                
                        the_idx = i + 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                        if (myid1 == the_id2) then
                            call MPI_SEND(By(1,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(By(Nx1,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                        end if
                
                    end if
                        
                    ! The Left condition
                    if (i == 1) then
                            
                        if (bcL == 1) then
                            the_idx = Nx_process
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(By(Nx,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(By(0,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                            end if
                        else if (bcL == 2) then
                            if (myid1 == the_id) then
                                By(0,0:Ny,k,d) = By(1,0:Ny,k,d)
                            end if
                        else if (bcL == 4) then
                            the_idx = Nx_process
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                            if (myid1 == the_id2) then
                                call MPI_SEND(By(Nx,0:Ny,k,d)*xa/xb,Ny + 1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(By(0,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if
                    else
                        the_idx = i - 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                        if (myid1 == the_id2) then
                            call MPI_SEND(By(Nx,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(By(0,0:Ny,k,d),Ny + 1,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                        end if
                    end if
                        
                end do
            end do
        end do
    end do
    
    end subroutine set_bc
    
    !*****************************************************************************************************
    
    subroutine calculate_pqs
    
    use com
    
    real(8) p1Gint3D(NumGLP,NumGLP,0:Nphi,NumEq),p2Gint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) q1Gint3D(NumGLP,NumGLP,0:Nphi,NumEq),q2Gint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) uhx1hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq),uhx2hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq)
    real(8) uhy1hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq),uhy2hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq)
    real(8) p1(Nx,Ny,0:Nphi,dimPk,NumEq),p2(Nx,Ny,0:Nphi,dimPk,NumEq)
    real(8) q1(Nx,Ny,0:Nphi,dimPk,NumEq),q2(Nx,Ny,0:Nphi,dimPk,NumEq)
    
    real(8),allocatable :: UR(:,:,:,:,:),UL(:,:,:,:,:),UU(:,:,:,:,:),UD(:,:,:,:,:)
    real(8),allocatable :: pR(:,:,:,:,:),pL(:,:,:,:,:),qU(:,:,:,:,:),qD(:,:,:,:,:)
    real(8),allocatable :: pU(:,:,:,:,:),pD(:,:,:,:,:),qR(:,:,:,:,:),qL(:,:,:,:,:)
    real(8),allocatable :: FR(:,:,:,:,:),FL(:,:,:,:,:),FU(:,:,:,:,:),FD(:,:,:,:,:)
    real(8),allocatable :: mu1R(:,:,:,:,:),mu1L(:,:,:,:,:),mu2U(:,:,:,:,:),mu2D(:,:,:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:),mu1hat(:,:,:,:,:),mu2hat(:,:,:,:,:)
    real(8),allocatable :: uhxhat(:,:,:,:,:),uhyhat(:,:,:,:,:),phxhat(:,:,:,:,:),qhyhat(:,:,:,:,:)
    real(8),allocatable :: ph(:,:,:,:,:),qh(:,:,:,:,:),uhpre(:,:,:,:,:)
    real(8),allocatable :: JzVertex(:,:,:)
    real(8),allocatable :: JzRL(:,:,:,:),JzUD(:,:,:,:)
    
    
    allocate(ph(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    allocate(qh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    allocate(uhpre(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    
    allocate(UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(UD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(pR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(pL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    
    allocate(pU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(pD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(qU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(qD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(FR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(FD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(mu1R(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu1L(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2U(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2D(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(JzRL(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(JzUD(Nx,0:Ny,0:Nphi,NumGLP))
    allocate(JzVertex(0:Nx,0:Ny,0:Nphi))
    
    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(phxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu1hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    
    
    uhpre = uh
    
    uh = uhs
    
    call set_bc
    call set_bc
    
    UR = 0
    UL = 0
    UU = 0
    UD = 0
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR(i,j,k,:,n) = UR(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR(:,d)
                        UL(i + 1,j,k,:,n) = UL(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 0,Ny
                    do i = 1,Nx
                        UU(i,j,k,:,n) = UU(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU(:,d)
                        UD(i,j + 1,k,:,n) = UD(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    URs = UR
    ULs = UL
    UUs = UU
    UDs = UD
    
    ! LDG scheme for ph and qh
    ! where ph approx uh_x, and qh approx uh_y
    ph = 0
    qh = 0
    p1 = 0
    p2 = 0
    q1 = 0
    q2 = 0    
    ! calculate the Volume integral
    do j = 1,Ny
        do i = 1,Nx
            
            uGint3D = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    do k = 0,Nphi
                        uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            do n = 1,NumEq
                do d = 1,dimPk1
                    do k = 0,Nphi
                        do j1 = 1,NumGLP
                            do i1 = 1,NumGLP
                                if (d > 1) then
                                    p1(i,j,k,d,n) = p1(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    p2(i,j,k,d,n) = p2(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    q1(i,j,k,d,n) = q1(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                    q2(i,j,k,d,n) = q2(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                    ph(i,j,k,d,n) = ph(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    qh(i,j,k,d,n) = qh(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    
    ! calculate the flux
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    UR1 = UL(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:)
                    uhx1hat(i,j,k,j1,:) = UR1
                    uhx2hat(i,j,k,j1,:) = UL1
                    uhxhat(i,j,k,j1,:) = 0.5*(UR1 + UL1)
                end do
            end do
        end do
    end do
    
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    UR1 = UD(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:)
                    uhy1hat(i,j,k,i1,:) = UR1
                    uhy2hat(i,j,k,i1,:) = UL1
                    uhyhat(i,j,k,i1,:) = 0.5*(UR1 + UL1)
                end do
            end do
        end do
    end do
    
    ! calculate the Surface integral
    do n = 1,NumEq
        do d = 1,dimPk1
            do j1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            p1(i,j,k,d,n) = p1(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhx1hat(i,j,k,j1,n)*phiGR(j1,d) - uhx1hat(i - 1,j,k,j1,n)*phiGL(j1,d))
                            p2(i,j,k,d,n) = p2(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhx2hat(i,j,k,j1,n)*phiGR(j1,d) - uhx2hat(i - 1,j,k,j1,n)*phiGL(j1,d))
                            ph(i,j,k,d,n) = ph(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhxhat(i,j,k,j1,n)*phiGR(j1,d) - uhxhat(i - 1,j,k,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk1
            do i1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            q1(i,j,k,d,n) = q1(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhy1hat(i,j,k,i1,n)*phiGU(i1,d) - uhy1hat(i,j - 1,k,i1,n)*phiGD(i1,d))
                            q2(i,j,k,d,n) = q2(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhy2hat(i,j,k,i1,n)*phiGU(i1,d) - uhy2hat(i,j - 1,k,i1,n)*phiGD(i1,d))
                            qh(i,j,k,d,n) = qh(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhyhat(i,j,k,i1,n)*phiGU(i1,d) - uhyhat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        p1(:,:,:,d,:) = p1(:,:,:,d,:)/mm(d)
        p2(:,:,:,d,:) = p2(:,:,:,d,:)/mm(d)
        q1(:,:,:,d,:) = q1(:,:,:,d,:)/mm(d)
        q2(:,:,:,d,:) = q2(:,:,:,d,:)/mm(d)
        ph(:,:,:,d,:) = ph(:,:,:,d,:)/mm(d)
        qh(:,:,:,d,:) = qh(:,:,:,d,:)/mm(d)
    end do
    
    uh = ph
    call set_bc
    call set_bc
    ph = uh
    
    uh = qh
    call set_bc
    call set_bc
    qh = uh
    
    phs = ph
    qhs = qh
    p1s = p1
    p2s = p2
    q1s = q1
    q2s = q2
    
    URUs = 0
    ULUs = 0
    URDs = 0
    ULDs = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do k = 0,Nphi
                do d = 1,dimPk
                    URUs(i,j,k,:) = URUs(i,j,k,:) + uh(i,j,k,d,:)*phiRU(d)
                    ULUs(i,j,k,:) = ULUs(i,j,k,:) + uh(i,j,k,d,:)*phiLU(d)
                    URDs(i,j,k,:) = URDs(i,j,k,:) + uh(i,j,k,d,:)*phiRD(d)
                    ULDs(i,j,k,:) = ULDs(i,j,k,:) + uh(i,j,k,d,:)*phiLD(d)
                end do
            end do
        end do
    end do
    
    uh = uhpre
    
    usGint3D = 0
    p1sGint3D = 0
    p2sGint3D = 0
    q1sGint3D = 0
    q2sGint3D = 0
    psGint3D = 0
    qsGint3D = 0
    do j = 1,Ny
        do i = 1,Nx
            
            do n = 1,NumEq
                do d = 1,dimPk
                    do k = 0,Nphi
                        usGint3D(:,:,i,j,k,n) = usGint3D(:,:,i,j,k,n) + uhs(i,j,k,d,n)*phiG(:,:,d)
                        p1sGint3D(:,:,i,j,k,n) = p1sGint3D(:,:,i,j,k,n) + p1s(i,j,k,d,n)*phiG(:,:,d)
                        p2sGint3D(:,:,i,j,k,n) = p2sGint3D(:,:,i,j,k,n) + p2s(i,j,k,d,n)*phiG(:,:,d)
                        q1sGint3D(:,:,i,j,k,n) = q1sGint3D(:,:,i,j,k,n) + q1s(i,j,k,d,n)*phiG(:,:,d)
                        q2sGint3D(:,:,i,j,k,n) = q2sGint3D(:,:,i,j,k,n) + q2s(i,j,k,d,n)*phiG(:,:,d)
                        psGint3D(:,:,i,j,k,n) = psGint3D(:,:,i,j,k,n) + phs(i,j,k,d,n)*phiG(:,:,d)
                        qsGint3D(:,:,i,j,k,n) = qsGint3D(:,:,i,j,k,n) + qhs(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
        end do
    end do
    
    end subroutine calculate_pqs
    
    !*****************************************************************************************************
    
    subroutine Lh
    
    use com
    
    real(8) Fx(NumGLP,NumGLP,0:Nphi,NumEq), Fy(NumGLP,NumGLP,0:Nphi,NumEq), Fz(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) rhoij,uij,vij,wij,Eij,B1ij,B2ij,B3ij,pij,Sij,Tij,Kij,rB1ij,rB2ij,rB3ij
    real(8) p1(Nx,Ny,0:Nphi,dimPk,NumEq),p2(Nx,Ny,0:Nphi,dimPk,NumEq)
    real(8) q1(Nx,Ny,0:Nphi,dimPk,NumEq),q2(Nx,Ny,0:Nphi,dimPk,NumEq)
    real(8) p1Gint3D(NumGLP,NumGLP,0:Nphi,NumEq),p2Gint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) q1Gint3D(NumGLP,NumGLP,0:Nphi,NumEq),q2Gint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) Hhat(NumGLP,NumGLP,0:Nphi,NumEq),uhij(NumEq)
    real(8) uhx1hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq),uhx2hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq)
    real(8) uhy1hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq),uhy2hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq)
    
    real(8),allocatable :: UR(:,:,:,:,:),UL(:,:,:,:,:),UU(:,:,:,:,:),UD(:,:,:,:,:)
    real(8),allocatable :: pR(:,:,:,:,:),pL(:,:,:,:,:),qU(:,:,:,:,:),qD(:,:,:,:,:)
    real(8),allocatable :: pU(:,:,:,:,:),pD(:,:,:,:,:),qR(:,:,:,:,:),qL(:,:,:,:,:)
    real(8),allocatable :: FR(:,:,:,:,:),FL(:,:,:,:,:),FU(:,:,:,:,:),FD(:,:,:,:,:)
    real(8),allocatable :: mu1R(:,:,:,:,:),mu1L(:,:,:,:,:),mu2U(:,:,:,:,:),mu2D(:,:,:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:),mu1hat(:,:,:,:,:),mu2hat(:,:,:,:,:)
    real(8),allocatable :: uhxhat(:,:,:,:,:),uhyhat(:,:,:,:,:),phxhat(:,:,:,:,:),qhyhat(:,:,:,:,:)
    real(8),allocatable :: ph(:,:,:,:,:),qh(:,:,:,:,:),uhpre(:,:,:,:,:)
    real(8),allocatable :: JzVertex(:,:,:)
    real(8),allocatable :: JzRL(:,:,:,:),JzUD(:,:,:,:)
    
    allocate(ph(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    allocate(qh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    allocate(uhpre(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq))
    
    allocate(UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(UD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(pR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(pL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    
    allocate(pU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(pD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(qU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(qD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(FR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(FD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(mu1R(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu1L(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2U(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2D(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    allocate(JzRL(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(JzUD(Nx,0:Ny,0:Nphi,NumGLP))
    allocate(JzVertex(0:Nx,0:Ny,0:Nphi))
    
    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(phxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu1hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    
    call set_bc
    call set_bc
    
    UR = 0
    UL = 0
    UU = 0
    UD = 0
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR(i,j,k,:,n) = UR(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR(:,d)
                        UL(i + 1,j,k,:,n) = UL(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 0,Ny
                    do i = 1,Nx
                        UU(i,j,k,:,n) = UU(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU(:,d)
                        UD(i,j + 1,k,:,n) = UD(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    ! LDG scheme for ph and qh
    ! where ph approx uh_x, and qh approx uh_y
    ph = 0
    qh = 0
    p1 = 0
    p2 = 0
    q1 = 0
    q2 = 0    
    ! calculate the Volume integral
    do j = 1,Ny
        do i = 1,Nx
            
            uGint3D = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    do k = 0,Nphi
                        uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            do n = 1,NumEq
                do d = 1,dimPk1
                    do k = 0,Nphi
                        do j1 = 1,NumGLP
                            do i1 = 1,NumGLP
                                if (d > 1) then
                                    p1(i,j,k,d,n) = p1(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    p2(i,j,k,d,n) = p2(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    q1(i,j,k,d,n) = q1(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                    q2(i,j,k,d,n) = q2(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                    ph(i,j,k,d,n) = ph(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phixG(i1,j1,d)
                                    qh(i,j,k,d,n) = qh(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*uGint3D(i1,j1,k,n)*phiyG(i1,j1,d)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    
    ! calculate the flux
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    UR1 = UL(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:)
                    uhx1hat(i,j,k,j1,:) = UR1
                    uhx2hat(i,j,k,j1,:) = UL1
                    uhxhat(i,j,k,j1,:) = 0.5*(UR1 + UL1)
                end do
            end do
        end do
    end do
    
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    UR1 = UD(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:)
                    uhy1hat(i,j,k,i1,:) = UR1
                    uhy2hat(i,j,k,i1,:) = UL1
                    uhyhat(i,j,k,i1,:) = 0.5*(UR1 + UL1)
                end do
            end do
        end do
    end do
    
    ! calculate the Surface integral
    do n = 1,NumEq
        do d = 1,dimPk1
            do j1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            p1(i,j,k,d,n) = p1(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhx1hat(i,j,k,j1,n)*phiGR(j1,d) - uhx1hat(i - 1,j,k,j1,n)*phiGL(j1,d))
                            p2(i,j,k,d,n) = p2(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhx2hat(i,j,k,j1,n)*phiGR(j1,d) - uhx2hat(i - 1,j,k,j1,n)*phiGL(j1,d))
                            ph(i,j,k,d,n) = ph(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(uhxhat(i,j,k,j1,n)*phiGR(j1,d) - uhxhat(i - 1,j,k,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk1
            do i1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            q1(i,j,k,d,n) = q1(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhy1hat(i,j,k,i1,n)*phiGU(i1,d) - uhy1hat(i,j - 1,k,i1,n)*phiGD(i1,d))
                            q2(i,j,k,d,n) = q2(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhy2hat(i,j,k,i1,n)*phiGU(i1,d) - uhy2hat(i,j - 1,k,i1,n)*phiGD(i1,d))
                            qh(i,j,k,d,n) = qh(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhyhat(i,j,k,i1,n)*phiGU(i1,d) - uhyhat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        p1(:,:,:,d,:) = p1(:,:,:,d,:)/mm(d)
        p2(:,:,:,d,:) = p2(:,:,:,d,:)/mm(d)
        q1(:,:,:,d,:) = q1(:,:,:,d,:)/mm(d)
        q2(:,:,:,d,:) = q2(:,:,:,d,:)/mm(d)
        ph(:,:,:,d,:) = ph(:,:,:,d,:)/mm(d)
        qh(:,:,:,d,:) = qh(:,:,:,d,:)/mm(d)
    end do
    
    uhpre = uh
    
    uh = ph
    call set_bc
    call set_bc
    ph = uh
    
    uh = qh
    call set_bc
    call set_bc
    qh = uh
    
    uh = uhpre
    
    ! DG scheme for uh        
    du = 0
    
    ! calculate the Volume integral
    do j = 1,Ny
        do i = 1,Nx
            
            uGint3D = 0
            pGint3D = 0
            qGint3D = 0
            p1Gint3D = 0
            p2Gint3D = 0
            q1Gint3D = 0
            q2Gint3D = 0
            do n = 1,NumEq
                do d = 1,dimPk
                    do k = 0,Nphi
                        uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                        p1Gint3D(:,:,k,n) = p1Gint3D(:,:,k,n) + p1(i,j,k,d,n)*phiG(:,:,d)
                        p2Gint3D(:,:,k,n) = p2Gint3D(:,:,k,n) + p2(i,j,k,d,n)*phiG(:,:,d)
                        q1Gint3D(:,:,k,n) = q1Gint3D(:,:,k,n) + q1(i,j,k,d,n)*phiG(:,:,d)
                        q2Gint3D(:,:,k,n) = q2Gint3D(:,:,k,n) + q2(i,j,k,d,n)*phiG(:,:,d)
                        pGint3D(:,:,k,n) = pGint3D(:,:,k,n) + ph(i,j,k,d,n)*phiG(:,:,d)
                        qGint3D(:,:,k,n) = qGint3D(:,:,k,n) + qh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            RHS = 0
            
            do k = 0,Nphi
                do j1 = 1,NumGLP
                    do i1 = 1,NumGLP
                        
                        uhij = uGint3D(i1,j1,k,:)
                        p1ij = p1Gint3D(i1,j1,k,:)
                        p2ij = p2Gint3D(i1,j1,k,:)
                        q1ij = q1Gint3D(i1,j1,k,:)
                        q2ij = q2Gint3D(i1,j1,k,:)
                        
                        uhsij = usGint3D(i1,j1,i,j,k,:)
                        p1sij = p1sGint3D(i1,j1,i,j,k,:)
                        p2sij = p2sGint3D(i1,j1,i,j,k,:)
                        q1sij = q1sGint3D(i1,j1,i,j,k,:)
                        q2sij = q2sGint3D(i1,j1,i,j,k,:)
                        
                        p1wij = p1ij - p1sij
                        p2wij = p2ij - p2sij
                        q1wij = q1ij - q1sij
                        q2wij = q2ij - q2sij
                        
                        call eigenvalueMm(SR,SL,uGint3D(i1,j1,k,1),uGint3D(i1,j1,k,2),uGint3D(i1,j1,k,3),uGint3D(i1,j1,k,4),uGint3D(i1,j1,k,5),uGint3D(i1,j1,k,6),uGint3D(i1,j1,k,7),uGint3D(i1,j1,k,8),1,0)
                        call eigenvalueMm(SU,SD,uGint3D(i1,j1,k,1),uGint3D(i1,j1,k,2),uGint3D(i1,j1,k,3),uGint3D(i1,j1,k,4),uGint3D(i1,j1,k,5),uGint3D(i1,j1,k,6),uGint3D(i1,j1,k,7),uGint3D(i1,j1,k,8),0,1)
                        
                        if (flux_type == 1) then
                            call calculate_H(Hij,0.5*(p1ij + p2ij),0.5*(q1ij + q2ij),uhij,0.5*(p1sij + p2sij),0.5*(q1sij + q2sij),uhsij,gamma)
                            call LF_Hamiltonian_2D
                        else if (flux_type == 2) then
                            call calculate_H(HRU,p1ij,q1ij,uhij,p1sij,q1sij,uhsij,gamma)
                            call calculate_H(HLU,p2ij,q1ij,uhij,p2sij,q1sij,uhsij,gamma)
                            call calculate_H(HRD,p1ij,q2ij,uhij,p1sij,q2sij,uhsij,gamma)
                            call calculate_H(HLD,p2ij,q2ij,uhij,p2sij,q2sij,uhsij,gamma)
                            call HLL_Hamiltonian_2D
                        end if
                        
                        Hhat(i1,j1,k,:) = Hhat1
                        
                        mu1(i1,j1,k,1) = 0
                        mu1(i1,j1,k,2) = pGint3D(i1,j1,k,2)*nu
                        mu1(i1,j1,k,3) = pGint3D(i1,j1,k,3)*nu
                        mu1(i1,j1,k,4) = pGint3D(i1,j1,k,4)*nu
                        mu1(i1,j1,k,5) = 0
                        mu1(i1,j1,k,6) = 0
                        mu1(i1,j1,k,7) = pGint3D(i1,j1,k,7)*eta - qGint3D(i1,j1,k,6)*eta
                        mu1(i1,j1,k,8) = 0
                        
                        mu2(i1,j1,k,1) = 0
                        mu2(i1,j1,k,2) = qGint3D(i1,j1,k,2)*nu
                        mu2(i1,j1,k,3) = qGint3D(i1,j1,k,3)*nu
                        mu2(i1,j1,k,4) = qGint3D(i1,j1,k,4)*nu
                        mu2(i1,j1,k,5) = 0
                        mu2(i1,j1,k,6) = qGint3D(i1,j1,k,6)*eta - pGint3D(i1,j1,k,7)*eta
                        mu2(i1,j1,k,7) = 0
                        mu2(i1,j1,k,8) = 0
                        
                        if (RHSCopen == 1) then
                            RHS(i1,j1,k,2) = RHS(i1,j1,k,2) + S2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k),tRK,eta,nu)
                            RHS(i1,j1,k,3) = RHS(i1,j1,k,3) + S3(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k),tRK,eta,nu)
                            RHS(i1,j1,k,5) = RHS(i1,j1,k,5) + S5(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k),tRK,gamma,eta,nu)
                        end if
                    end do
                end do
            end do
            
            do n = 1,NumEq
                do d = 1,dimPk1
                    do k = 0,Nphi
                        do j1 = 1,NumGLP
                            do i1 = 1,NumGLP
                                du(i,j,k,d,n) = du(i,j,k,d,n) - 0.25d0*weight(i1)*weight(j1)*(mu1(i1,j1,k,n)*phixG(i1,j1,d) + mu2(i1,j1,k,n)*phiyG(i1,j1,d) + Hhat(i1,j1,k,n)*phiG(i1,j1,d))
                                du(i,j,k,d,n) = du(i,j,k,d,n) + 0.25d0*weight(i1)*weight(j1)*RHS(i1,j1,k,n)*phiG(i1,j1,d)
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    
    ! calculate the Numerical flux
    pR = 0
    pL = 0
    pU = 0
    pD = 0
    
    qR = 0
    qL = 0
    qU = 0
    qD = 0
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        pR(i,j,k,:,n) = pR(i,j,k,:,n) + ph(i,j,k,d,n)*phiGR(:,d)
                        pL(i + 1,j,k,:,n) = pL(i + 1,j,k,:,n) + ph(i + 1,j,k,d,n)*phiGL(:,d)
                        qR(i,j,k,:,n) = qR(i,j,k,:,n) + qh(i,j,k,d,n)*phiGR(:,d)
                        qL(i + 1,j,k,:,n) = qL(i + 1,j,k,:,n) + qh(i + 1,j,k,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 0,Ny
                    do i = 1,Nx
                        qU(i,j,k,:,n) = qU(i,j,k,:,n) + qh(i,j,k,d,n)*phiGU(:,d)
                        qD(i,j + 1,k,:,n) = qD(i,j + 1,k,:,n) + qh(i,j + 1,k,d,n)*phiGD(:,d)
                        pU(i,j,k,:,n) = pU(i,j,k,:,n) + ph(i,j,k,d,n)*phiGU(:,d)
                        pD(i,j + 1,k,:,n) = pD(i,j + 1,k,:,n) + ph(i,j + 1,k,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do
    
    ! The x-flux
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    
                    uij = uR(i,j,k,j1,2)
                    vij = uR(i,j,k,j1,3)
                    B1ij = uR(i,j,k,j1,6)
                    B2ij = uR(i,j,k,j1,7)
                    
                    FR(i,j,k,j1,6) = 0
                    FR(i,j,k,j1,7) = uij*B2ij - vij*B1ij
                    
                    mu1R(i,j,k,j1,1) = 0
                    mu1R(i,j,k,j1,2) = pR(i,j,k,j1,2)*nu
                    mu1R(i,j,k,j1,3) = pR(i,j,k,j1,3)*nu
                    mu1R(i,j,k,j1,4) = pR(i,j,k,j1,4)*nu
                    mu1R(i,j,k,j1,5) = 0
                    mu1R(i,j,k,j1,6) = 0
                    mu1R(i,j,k,j1,7) = pR(i,j,k,j1,7)*eta - qR(i,j,k,j1,6)*eta
                    mu1R(i,j,k,j1,8) = 0
                end do
            end do
        end do
    end do
    
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 1,Nx1
                    
                    uij = uL(i,j,k,j1,2)
                    vij = uL(i,j,k,j1,3)
                    B1ij = uL(i,j,k,j1,6)
                    B2ij = uL(i,j,k,j1,7)
                    
                    FL(i,j,k,j1,6) = 0
                    FL(i,j,k,j1,7) = uij*B2ij - vij*B1ij
                    
                    mu1L(i,j,k,j1,1) = 0
                    mu1L(i,j,k,j1,2) = pL(i,j,k,j1,2)*nu
                    mu1L(i,j,k,j1,3) = pL(i,j,k,j1,3)*nu
                    mu1L(i,j,k,j1,4) = pL(i,j,k,j1,4)*nu
                    mu1L(i,j,k,j1,5) = 0
                    mu1L(i,j,k,j1,6) = 0
                    mu1L(i,j,k,j1,7) = pL(i,j,k,j1,7)*eta - qL(i,j,k,j1,6)*eta
                    mu1L(i,j,k,j1,8) = 0
                end do
            end do
        end do
    end do
    
    ! The y-Flux
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    
                    uij = UU(i,j,k,i1,2)
                    vij = UU(i,j,k,i1,3)
                    B1ij = UU(i,j,k,i1,6)
                    B2ij = UU(i,j,k,i1,7)
                    
                    FU(i,j,k,i1,6) = vij*B1ij - uij*B2ij
                    FU(i,j,k,i1,7) = 0
                    
                    mu2U(i,j,k,i1,1) = 0
                    mu2U(i,j,k,i1,2) = qU(i,j,k,i1,2)*nu
                    mu2U(i,j,k,i1,3) = qU(i,j,k,i1,3)*nu
                    mu2U(i,j,k,i1,4) = qU(i,j,k,i1,4)*nu
                    mu2U(i,j,k,i1,5) = 0
                    mu2U(i,j,k,i1,6) = qU(i,j,k,i1,6)*eta - pU(i,j,k,i1,7)*eta
                    mu2U(i,j,k,i1,7) = 0
                    mu2U(i,j,k,i1,8) = 0
                end do
            end do
        end do
    end do
    
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny1
                do i = 1,Nx
                    
                    uij = UD(i,j,k,i1,2)
                    vij = UD(i,j,k,i1,3)
                    B1ij = UD(i,j,k,i1,6)
                    B2ij = UD(i,j,k,i1,7)
                    
                    FD(i,j,k,i1,6) = vij*B1ij - uij*B2ij
                    FD(i,j,k,i1,7) = 0
                    
                    mu2D(i,j,k,i1,1) = 0
                    mu2D(i,j,k,i1,2) = qD(i,j,k,i1,2)*nu
                    mu2D(i,j,k,i1,3) = qD(i,j,k,i1,3)*nu
                    mu2D(i,j,k,i1,4) = qD(i,j,k,i1,4)*nu
                    mu2D(i,j,k,i1,5) = 0
                    mu2D(i,j,k,i1,6) = qD(i,j,k,i1,6)*eta - pD(i,j,k,i1,7)*eta
                    mu2D(i,j,k,i1,7) = 0
                    mu2D(i,j,k,i1,8) = 0
                end do
            end do
        end do
    end do
    
    ! calculate Fx hat
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    call eigenvalueMm(SRmax,SRmin,UR(i,j,k,j1,1),UR(i,j,k,j1,2),UR(i,j,k,j1,3),UR(i,j,k,j1,4),UR(i,j,k,j1,5),UR(i,j,k,j1,6),UR(i,j,k,j1,7),UR(i,j,k,j1,8),1,0)
                    call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,k,j1,1),UL(i + 1,j,k,j1,2),UL(i + 1,j,k,j1,3),UL(i + 1,j,k,j1,4),UL(i + 1,j,k,j1,5),UL(i + 1,j,k,j1,6),UL(i + 1,j,k,j1,7),UL(i + 1,j,k,j1,8),1,0)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FL(i + 1,j,k,j1,:)
                    FL1 = FR(i,j,k,j1,:)
                    UR1 = UL(i + 1,j,k,j1,:) - ULs(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:) - URs(i,j,k,j1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    end if
                    Fxhat(i,j,k,j1,:) = Fhat1
                    mu1hat(i,j,k,j1,:) = 0.5*(mu1R(i,j,k,j1,:) + mu1L(i + 1,j,k,j1,:))
                end do
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    call eigenvalueMm(SRmax,SRmin,UU(i,j,k,i1,1),UU(i,j,k,i1,2),UU(i,j,k,i1,3),UU(i,j,k,i1,4),UU(i,j,k,i1,5),UU(i,j,k,i1,6),UU(i,j,k,i1,7),UU(i,j,k,i1,8),0,1)
                    call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,k,i1,1),UD(i,j + 1,k,i1,2),UD(i,j + 1,k,i1,3),UD(i,j + 1,k,i1,4),UD(i,j + 1,k,i1,5),UD(i,j + 1,k,i1,6),UD(i,j + 1,k,i1,7),UD(i,j + 1,k,i1,8),0,1)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FD(i,j + 1,k,i1,:)
                    FL1 = FU(i,j,k,i1,:)
                    UR1 = UD(i,j + 1,k,i1,:) - UDs(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:) - UUs(i,j,k,i1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    end if
                    Fyhat(i,j,k,i1,:) = Fhat1
                    mu2hat(i,j,k,i1,:) = 0.5*(mu2U(i,j,k,i1,:) + mu2D(i,j + 1,k,i1,:))
                end do
            end do
        end do
    end do
    
    ! calculate the Surface integral
    do n = 1,NumEq
        do d = 1,dimPk1
            do j1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            du(i,j,k,d,n) = du(i,j,k,d,n) + (0.5d0/hx)*weight(j1)*(mu1hat(i,j,k,j1,n)*phiGR(j1,d) - mu1hat(i - 1,j,k,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do n = 1,NumEq
        do d = 1,dimPk1
            do i1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            du(i,j,k,d,n) = du(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(mu2hat(i,j,k,i1,n)*phiGU(i1,d) - mu2hat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        du(:,:,:,d,:) = du(:,:,:,d,:)/mm(d)
    end do
    
    ! Scheme of Bx and By
    EzRL = -Fxhat(:,:,:,:,7)
    EzUD = Fyhat(:,:,:,:,6)
    JzRL = mu1hat(:,:,:,:,7)
    JzUD = -mu2hat(:,:,:,:,6)
    
    ! calculate each Uh and Ph and Qh at vertex
    URU = 0
    ULU = 0
    URD = 0
    ULD = 0
    PRU = 0
    PLU = 0
    PRD = 0
    PLD = 0
    QRU = 0
    QLU = 0
    QRD = 0
    QLD = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do k = 0,Nphi
                do d = 1,dimPk
                    URU(i,j,k,:) = URU(i,j,k,:) + uh(i,j,k,d,:)*phiRU(d)
                    ULU(i,j,k,:) = ULU(i,j,k,:) + uh(i,j,k,d,:)*phiLU(d)
                    URD(i,j,k,:) = URD(i,j,k,:) + uh(i,j,k,d,:)*phiRD(d)
                    ULD(i,j,k,:) = ULD(i,j,k,:) + uh(i,j,k,d,:)*phiLD(d)
                    
                    PRU(i,j,k,:) = PRU(i,j,k,:) + ph(i,j,k,d,:)*phiRU(d)
                    PLU(i,j,k,:) = PLU(i,j,k,:) + ph(i,j,k,d,:)*phiLU(d)
                    PRD(i,j,k,:) = PRD(i,j,k,:) + ph(i,j,k,d,:)*phiRD(d)
                    PLD(i,j,k,:) = PLD(i,j,k,:) + ph(i,j,k,d,:)*phiLD(d)
                    
                    QRU(i,j,k,:) = QRU(i,j,k,:) + qh(i,j,k,d,:)*phiRU(d)
                    QLU(i,j,k,:) = QLU(i,j,k,:) + qh(i,j,k,d,:)*phiLU(d)
                    QRD(i,j,k,:) = QRD(i,j,k,:) + qh(i,j,k,d,:)*phiRD(d)
                    QLD(i,j,k,:) = QLD(i,j,k,:) + qh(i,j,k,d,:)*phiLD(d)
                end do
            end do
        end do
    end do
    
    ! calculate Ez and Jz at vertex
    do i = 0,Nx
        do j = 0,Ny
            do k = 0,Nphi
                URU1 = ULD(i + 1,j + 1,k,:)
                ULU1 = URD(i,j + 1,k,:)
                URD1 = ULU(i + 1,j,k,:)
                ULD1 = URU(i,j,k,:)
                
                URUs1 = ULDs(i + 1,j + 1,k,:)
                ULUs1 = URDs(i,j + 1,k,:)
                URDs1 = ULUs(i + 1,j,k,:)
                ULDs1 = URUs(i,j,k,:)
                
                PRU1 = PLD(i + 1,j + 1,k,:)
                PLU1 = PRD(i,j + 1,k,:)
                PRD1 = PLU(i + 1,j,k,:)
                PLD1 = PRU(i,j,k,:)
                
                QRU1 = QLD(i + 1,j + 1,k,:)
                QLU1 = QRD(i,j + 1,k,:)
                QRD1 = QLU(i + 1,j,k,:)
                QLD1 = QRU(i,j,k,:)
                
                EEzRU = PRU1(7) - QRU1(6)  
                EEzLU = PLU1(7) - QLU1(6) 
                EEzRD = PRD1(7) - QRD1(6) 
                EEzLD = PLD1(7) - QLD1(6) 
            
                if (flux_type == 1) then
                    call LF_Flux_2D
                else if (flux_type == 2) then
                    call HLL_Flux_2D
                end if
            
                EzVertex(i,j,k) = Ezhat
                JzVertex(i,j,k) = 0.25*eta*(EEzRU + EEzLU + EEzRD + EEzLD)
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
                        dBx(i,j,k,d) = dBx(i,j,k,d) + 0.5*weight(j1)*(EzRL(i,j,k,j1) + JzRL(i,j,k,j1))*EzyG(j1,d)
                    end do
                    dBx(i,j,k,d) = (dBx(i,j,k,d) - (EzVertex(i,j,k) + JzVertex(i,j,k))*EzU(d)/hy + (EzVertex(i,j - 1,k) + JzVertex(i,j - 1,k))*EzD(d)/hy)/mmE(d)
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
                        dBy(i,j,k,d) = dBy(i,j,k,d) - 0.5*weight(i1)*(EzUD(i,j,k,i1) + JzUD(i,j,k,i1))*EzxG(i1,d)
                    end do
                    dBy(i,j,k,d) = (dBy(i,j,k,d) + (EzVertex(i,j,k) + JzVertex(i,j,k))*EzR(d)/hx - (EzVertex(i - 1,j,k) + JzVertex(i - 1,j,k))*EzL(d)/hx)/mmE(d)
                end do
            end do
        end do
    end do
    
    !uh(1:Nx,1:Ny,:,:,:) = p1 - p1s
    !uh = du
    
    end subroutine Lh
    
    !*****************************************************************************************************
    
    subroutine LF_Flux
    
    use com
    
    Fhat1 = 0.5d0*(FR1 + FL1 - max(abs(SR),abs(SL))*(UR1 - UL1))
    
    end subroutine LF_Flux
    
    !*****************************************************************************************************
    
    subroutine div_free
    
    use com
    
    real BxR(kk + 1),BxL(kk + 1),ByU(kk + 1),ByD(kk + 1)
    real Bxint(dimPk),Byint(dimPk)
    real a0R,a1R,a2R,a0L,a1L,a2L
    real b0U,b1U,b2U,b0D,b1D,b2D
    real c00,c10,c01,c20,c11,c02
    real a00,a10,a01,a20,a11,a02,a30,a21,a12
    real b00,b10,b01,b20,b11,b02,b21,b12,b03
    real c5,c7,c8
    real a00new,a10new,a01new,a20new,a11new,a02new
    real b00new,b10new,b01new,b20new,b11new,b02new
    real rxy,ryx,ri,S1,S2
    
    !uh(:,:,:,4:dimPk,:) = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
            
                a00 = uh(i,j,k,1,1)
                a10 = uh(i,j,k,2,1)
                a01 = uh(i,j,k,3,1)
                a20 = uh(i,j,k,4,1)
                a11 = uh(i,j,k,5,1)
                a02 = uh(i,j,k,6,1)
                
                b00 = uh(i,j,k,1,2)
                b10 = uh(i,j,k,2,2)
                b01 = uh(i,j,k,3,2)
                b20 = uh(i,j,k,4,2)
                b11 = uh(i,j,k,5,2)
                b02 = uh(i,j,k,6,2)
    
                ! The reconstruction of B = (Bx,By) from the interface
                a00new = a00
                a01new = a01
                a02new = a02
                
                b00new = b00
                b10new = b10
                b20new = b20
                
                c5 = (a10*hx - b01*hy)/(hx**2 + hy**2)
                c7 = (2*a20*hx - 5*b11*hy)/(2*hx**2 + 10*hy**2) 
                c8 = (5*a11*hx - 2*b02*hy)/(10*hx**2 + 2*hy**2)
                
                a10new = hx*c5
                a20new = hx*c7
                a11new = 2*hx*c8
                
                b01new = -hy*c5
                b11new = -2*hy*c7
                b02new = -hy*c8
                
                !uh(i,j,k,1,1) = a00new
                !uh(i,j,k,2,1) = a10new
                !uh(i,j,k,3,1) = a01new
                !uh(i,j,k,4,1) = a20new
                !uh(i,j,k,5,1) = a11new
                !uh(i,j,k,6,1) = a02new
            
                !uh(i,j,k,1,2) = b00new
                !uh(i,j,k,2,2) = b10new
                !uh(i,j,k,3,2) = b01new
                !uh(i,j,k,4,2) = b20new
                !uh(i,j,k,5,2) = b11new
                !uh(i,j,k,6,2) = b02new
            
            end do
        end do
    end do
    
    end subroutine div_free
    
    !*****************************************************************************************************
    
    subroutine calculate_totaldiv
    
    use com
    real(8) uGdiv(NumGLP,NumGLP),totaldiv1
    
    totaldiv = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                uGdiv = 0
                do d = 1,dimPk
                    uGdiv = uGdiv + uh(i,j,k,d,6)*phixG(:,:,d) + uh(i,j,k,d,7)*phiyG(:,:,d)
                end do
            
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        totaldiv = totaldiv + 0.25*weight(i1)*weight(j1)*abs(uGdiv(i1,j1))
                        !if (abs(uGdiv(i1,j1)) > 1e-15) then
                        !    print *,uGdiv(i1,j1)
                        !end if
                    end do
                end do
            end do
        end do
    end do
    
    do the_id = 2,N_process
        
        if (myid1 == the_id) then
            call MPI_SEND(totaldiv,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if
        
        if (myid1 == 1) then
            call MPI_RECV(totaldiv1,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            totaldiv = totaldiv + totaldiv1
        end if
        
    end do
    
    if (myid1 == 1) then
        totaldiv = totaldiv/(Nx0*Ny0*Nphi1)
    end if
    
    end subroutine calculate_totaldiv
    
    !*****************************************************************************************************
    
    subroutine div_free_Balsara
    
    use com
    
    real(8) BxR(kk + 1),BxL(kk + 1),ByU(kk + 1),ByD(kk + 1)
    real(8) a0R,a1R,a2R,a0L,a1L,a2L
    real(8) b0U,b1U,b2U,b0D,b1D,b2D
    real(8) a00,a10,a01,a20,a11,a02,a30,a21,a12,a03,a40,a13
    real(8) b00,b10,b01,b20,b11,b02,b30,b21,b12,b03,b31,b04
    real(8) rxy,ryx,S1,S2,omega
    real(8) Bxbar1(0:Nx,0:Ny1,0:Nphi),Bybar1(0:Nx1,0:Ny,0:Nphi)
    
    omega1 = uh(1:Nx,1:Ny,:,2,7) - uh(1:Nx,1:Ny,:,3,6)
    uh(:,:,:,:,6:7) = 0
    
    Bxbar1 = Bx(:,:,:,1)
    Bybar1 = By(:,:,:,1)
    
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                if ( abs(-(Bx(i,j,k,1) - Bx(i - 1,j,k,1))*(hy/hx) + Bybar1(i,j - 1,k) - Bybar1(i,j,k)) >= 1e-15 ) then
                    Bybar1(i,j,k) = -(Bx(i,j,k,1) - Bx(i - 1,j,k,1))*(hy/hx) + Bybar1(i,j - 1,k)
                end if
                if ( abs(-(By(i,j,k,1) - By(i,j - 1,k,1))*(hx/hy) + Bxbar1(i - 1,j,k) - Bxbar1(i,j,k)) >= 1e-15 ) then
                    Bxbar1(i,j,k) = -(By(i,j,k,1) - By(i,j - 1,k,1))*(hx/hy) + Bxbar1(i - 1,j,k)
                end if
            end do
        end do
    end do
    
    Bx(:,:,:,1) = 0.5d0*(Bxbar1 + Bx(:,:,:,1))
    By(:,:,:,1) = 0.5d0*(Bybar1 + By(:,:,:,1))
    
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
                a3R = BxR(4)
    
                a0L = BxL(1)
                a1L = BxL(2)
                a2L = BxL(3)
                a3L = BxL(4)
    
                b0U = ByU(1)
                b1U = ByU(2)
                b2U = ByU(3)
                b3U = ByU(4)
    
                b0D = ByD(1)
                b1D = ByD(2)
                b2D = ByD(3)
                b3D = ByD(4)
            
                rxy = hx/hy
                ryx = hy/hx
                S1 = (1d0/2d0)*(a1R + a1L)
                S2 = (1d0/2d0)*(b1U + b1D)
                omega = omega1(i,j,k)
    
                ! The reconstruction of B = (Bx,By) from the interface
                a03 = (1d0/2d0)*(a3R + a3L)
                b30 = (1d0/2d0)*(b3U + b3D)
                a13 = (1d0/2d0)*(a3R - a3L)
                b31 = (1d0/2d0)*(b3U - b3D)
                a02 = (1d0/2d0)*(a2R + a2L)
                b20 = (1d0/2d0)*(b2U + b2D)
                a12 = (1d0/2d0)*(a2R - a2L)
                b21 = (1d0/2d0)*(b2U - b2D)
                a11 = (1d0/2d0)*(a1R - a1L)
                b11 = (1d0/2d0)*(b1U - b1D)
                a40 = -(1d0/4d0)*rxy*b31
                b04 = -(1d0/4d0)*ryx*a13
                a30 = -(1d0/3d0)*rxy*b21
                b03 = -(1d0/3d0)*ryx*a12
                a10 = (1d0/2d0)*(a0R - a0L) - (2d0/5d0)*a30
                b01 = (1d0/2d0)*(b0U - b0D) - (2d0/5d0)*b03
                a20 = -(1d0/2d0)*rxy*(b11 - 3d0/5d0*b31) + (6d0/7d0)*a40
                b02 = -(1d0/2d0)*ryx*(a11 - 3d0/5d0*a13) + (6d0/7d0)*b04
                a00 = (1d0/2d0)*(a0R + a0L) - (2d0/3d0)*a20 - (8d0/35d0)*a40
                b00 = (1d0/2d0)*(b0U + b0D) - (2d0/3d0)*b02 - (8d0/35d0)*b04
                a01 = (S1*ryx + S2 - omega)/(1 + ryx)
                b10 = a01 + omega
                a21 = (3d0/2d0)*(S1 - a01)
                b12 = (3d0/2d0)*(S2 - b10)
                
                !a00 = (1d0/2d0)*(a0R + a0L) + (1d0/6d0)*rxy*(b1U - b1D)
                !a10 = (1d0/2d0)*(a0R - a0L) + (1d0/15d0)*rxy*(b2U - b2D)
                !a01 = (1d0/2d0)*(a1R + a1L)
                !a20 = -(1d0/4d0)*rxy*(b1U - b1D)
                !a11 = (1d0/2d0)*(a1R - a1L)
                !a02 = (1d0/2d0)*(a2R + a2L)
                !a30 = -(1d0/6d0)*rxy*(b2U - b2D)
                !a12 = (1d0/2d0)*(a2R - a2L)
            
                !b00 = (1d0/2d0)*(b0U + b0D) + (1d0/6d0)*ryx*(a1R - a1L)
                !b01 = (1d0/2d0)*(b0U - b0D) + (1d0/15d0)*ryx*(a2R - a2L)
                !b10 = (1d0/2d0)*(b1U + b1D)
                !b02 = -(1d0/4d0)*ryx*(a1R - a1L)
                !b11 = (1d0/2d0)*(b1U - b1D)
                !b20 = (1d0/2d0)*(b2U + b2D)
                !b03 = -(1d0/6d0)*ryx*(a2R - a2L)
                !b21 = (1d0/2d0)*(b2U - b2D)
            
                uh(i,j,k,1,6) = a00
                uh(i,j,k,2,6) = a10
                uh(i,j,k,3,6) = a01
                uh(i,j,k,4,6) = a20
                uh(i,j,k,5,6) = a11
                uh(i,j,k,6,6) = a02
                uh(i,j,k,7,6) = a30
                uh(i,j,k,8,6) = a21
                uh(i,j,k,9,6) = a12
                uh(i,j,k,10,6) = a03
                uh(i,j,k,11,6) = a40
                uh(i,j,k,14,6) = a13
            
                uh(i,j,k,1,7) = b00
                uh(i,j,k,2,7) = b10
                uh(i,j,k,3,7) = b01
                uh(i,j,k,4,7) = b20
                uh(i,j,k,5,7) = b11
                uh(i,j,k,6,7) = b02
                uh(i,j,k,7,7) = b30
                uh(i,j,k,8,7) = b21
                uh(i,j,k,9,7) = b12
                uh(i,j,k,10,7) = b03
                uh(i,j,k,12,7) = b31
                uh(i,j,k,15,7) = b04
            
            end do
        end do
    end do
    
    end subroutine div_free_Balsara
    
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
    
    call eigenvalueMm(alphaxRU,betaxRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),1,0)
    call eigenvalueMm(alphayRU,betayRU,URU1(1),URU1(2),URU1(3),URU1(4),URU1(5),URU1(6),URU1(7),URU1(8),0,1)
    call eigenvalueMm(alphaxLU,betaxLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),1,0)
    call eigenvalueMm(alphayLU,betayLU,ULU1(1),ULU1(2),ULU1(3),ULU1(4),ULU1(5),ULU1(6),ULU1(7),ULU1(8),0,1)
    call eigenvalueMm(alphaxRD,betaxRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),1,0)
    call eigenvalueMm(alphayRD,betayRD,URD1(1),URD1(2),URD1(3),URD1(4),URD1(5),URD1(6),URD1(7),URD1(8),0,1)
    call eigenvalueMm(alphaxLD,betaxLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),1,0)
    call eigenvalueMm(alphayLD,betayLD,ULD1(1),ULD1(2),ULD1(3),ULD1(4),ULD1(5),ULD1(6),ULD1(7),ULD1(8),0,1)
    alphax2D = max(alphaxRU,alphaxLU,alphaxRD,alphaxLD)
    alphay2D = max(alphayRU,alphayLU,alphayRD,alphayLD)
    
    call calculate_Ez(EzRU,URU1(2),URU1(3),URU1(6),URU1(7))
    call calculate_Ez(EzLU,ULU1(2),ULU1(3),ULU1(6),ULU1(7))
    call calculate_Ez(EzRD,URD1(2),URD1(3),URD1(6),URD1(7))
    call calculate_Ez(EzLD,ULD1(2),ULD1(3),ULD1(6),ULD1(7))
    
    B1RU = URU1(6)
    B1LU = ULU1(6)
    B1RD = URD1(6)
    B1LD = ULD1(6)
    
    B2RU = URU1(7)
    B2LU = ULU1(7)
    B2RD = URD1(7)
    B2LD = ULD1(7)
    
    Ezhat = 0.25*(EzRU + EzLU + EzRD + EzLD) - 0.5*alphay2D*(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)) + 0.5*alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    
    end subroutine LF_Flux_2D
    
    !*****************************************************************************************************
    
    subroutine calculate_Ez(Ez,u1,u2,B1,B2)
    
    real(8) Ez,u1,u2,B1,B2
    
    Ez = u2*B1 - u1*B2
    
    end subroutine calculate_Ez
    
    !*****************************************************************************************************
    
    function S2(r,z,phi,t,eta,nu)
    
    real r,z,phi,t,eta,nu
    real S2
    parameter(pi = 4*atan(1d0))
	
    S2 = -nu*(1d0/(2*pi))*(4 - (r - t)**2 - (z - t)**2)*(z - t)*exp(0.5*(1 - (r - t)**2 - (z - t)**2))
    
	end
	!*************************************************
	function S3(r,z,phi,t,eta,nu)
    
    real r,z,phi,t,eta,nu
    real S3
    parameter(pi = 4*atan(1d0))
	
    S3 = nu*(1d0/(2*pi))*(4 - (r - t)**2 - (z - t)**2)*(r - t)*exp(0.5*(1 - (r - t)**2 - (z - t)**2))
    
    end
	!*************************************************
	function S5(r,z,phi,t,gamma,eta,nu)
    
    real r,z,phi,t,eta,nu
    real S5
	
    S5 = 0
    
    end
    
    !*****************************************************************************************************
    
    subroutine eigenvalueMm(Amax,Amin,rho,u,v,w,press,B1,B2,B3,n1,n2)
    
    use com
    
    real(8) u,v,w,p,c,BP,Bn,un,cf,n3,press
    
    n3 = 0
    
    BP = B1**2 + B2**2 + B3**2
    Bn = B1*n1 + B2*n2 + B3*n3
    un = u*n1 + v*n2 + w*n3
    
    p = press
    
    c = sqrt(abs(gamma*p/rho))
    
    cf = sqrt(abs( 0.5d0*(c**2 + BP/rho + sqrt((c**2 + BP/rho)**2 - 4*c**2*Bn**2/rho) ) ))
    
    Amax = un + cf
    Amin = un - cf
    
    end subroutine eigenvalueMm
    
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
    
    BxRU = URU1(6) - URUs1(6)
    BxLU = ULU1(6) - ULUs1(6)
    BxRD = URD1(6) - URDs1(6)
    BxLD = ULD1(6) - ULDs1(6)
    
    ByRU = URU1(7) - URUs1(7)
    ByLU = ULU1(7) - ULUs1(7)
    ByRD = URD1(7) - URDs1(7)
    ByLD = ULD1(7) - ULDs1(7)
    
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
    u = Uh(2)
    v = Uh(3)
    w = Uh(4)
    p = Uh(5)
    B1 = Uh(6)
    B2 = Uh(7)
    B3 = Uh(8)

    E = p/gamma1 + 0.5*rho*(u**2 + v**2 + w**2) + 0.5*(B1**2 + B2**2 + B3**2)
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
    
    subroutine evenex_y(a,b)
    
    real a(10),b(10)
    
    a(1) = b(1)
    
    a(2) = b(2)
    a(3) = -b(3)
    
    a(4) = b(4)
    a(5) = -b(5)
    a(6) = b(6)
    
    a(7) = b(7)
    a(8) = -b(8)
    a(9) = b(9)
    a(10) = -b(10)
    
    a(11) = b(11)
    a(12) = -b(12)
    a(13) = b(13)
    a(14) = -b(14)
    a(15) = b(15)
    
    end subroutine evenex_y
    
    !*****************************************************************************************************
    
    subroutine oddex_y(a,b)
    
    real a(10),b(10)
    
    a(1) = -b(1)
    
    a(2) = -b(2)
    a(3) = b(3)
    
    a(4) = -b(4)
    a(5) = b(5)
    a(6) = -b(6)
    
    a(7) = -b(7)
    a(8) = b(8)
    a(9) = -b(9)
    a(10) = b(10)
    
    a(11) = -b(11)
    a(12) = b(12)
    a(13) = -b(13)
    a(14) = b(14)
    a(15) = -b(15)
    
    end subroutine oddex_y
    
    !*****************************************************************************************************
    
    subroutine evenex_x(a,b)
    
    real a(10),b(10)
    
    a(1) = b(1)
    
    a(2) = -b(2)
    a(3) = b(3)
    
    a(4) = b(4)
    a(5) = -b(5)
    a(6) = b(6)
    
    a(7) = -b(7)
    a(8) = b(8)
    a(9) = -b(9)
    a(10) = b(10)
    
    a(11) = -b(11)
    a(12) = b(12)
    a(13) = -b(13)
    a(14) = b(14)
    a(15) = -b(15)
    
    end subroutine evenex_x
    
    !*****************************************************************************************************
    
    subroutine oddex_x(a,b)
    
    real a(10),b(10)
    
    a(1) = -b(1)
    
    a(2) = b(2)
    a(3) = -b(3)
    
    a(4) = -b(4)
    a(5) = b(5)
    a(6) = -b(6)
    
    a(7) = b(7)
    a(8) = -b(8)
    a(9) = b(9)
    a(10) = -b(10)
    
    a(11) = b(11)
    a(12) = -b(12)
    a(13) = b(13)
    a(14) = -b(14)
    a(15) = b(15)
    
    end subroutine oddex_x
    
    !*****************************************************************************************************
    
    subroutine save_solution2
    
    use com
    
    real uhsave(NumEq)
    integer the_idx1,the_idy1
    
    !uh = du
    uh = uhs
    
    if (myid1 == 1) then
        
        open(unit = 1,file = 'Q1s.txt')
        open(unit = 2,file = 'Q2s.txt')
        open(unit = 3,file = 'Q3s.txt')
        open(unit = 4,file = 'Q4s.txt')
        open(unit = 5,file = 'Q5s.txt')
        open(unit = 6,file = 'Q6s.txt')
        open(unit = 7,file = 'Q7s.txt')
        open(unit = 8,file = 'Q8s.txt')
        
    end if
    
    !do d = 1,dimPk
        do j = 1,Ny0
            do i = 1,Nx0
                    
                the_idx1 = mod(i,Nx)
                if (the_idx1 == 0) then
                    the_idx1 = Nx
                end if
                the_idx = (i - the_idx1)/Nx + 1
            
                the_idy1 = mod(j,Ny)
                if (the_idy1 == 0) then
                    the_idy1 = Ny
                end if
                the_idy = (j - the_idy1)/Ny + 1
            
                the_id = the_idx + Nx_process*(the_idy - 1)
                    
                if (the_id /= 1) then
                    if (myid1 == the_id) then
                        call MPI_SEND(uh(the_idx1,the_idy1,0,1,:),NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
                    end if
                    if (myid1 == 1) then
                        call MPI_RECV(uhsave,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
                    end if
                else if (the_id == 1) then
                    if (myid1 == 1) then
                        uhsave = uh(the_idx1,the_idy1,0,1,:)
                    end if
                end if
            
                if (myid1 == 1) then
                    write(1,*) uhsave(1)
                    write(2,*) uhsave(2)
                    write(3,*) uhsave(3)
                    write(4,*) uhsave(4)
                    write(5,*) uhsave(5)
                    write(6,*) uhsave(6)
                    write(7,*) uhsave(7)
                    write(8,*) uhsave(8)
                end if
                    
            end do
        end do
    !end do
    
    if (myid1 == 1) then
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
    end if
    
    end subroutine save_solution2
    
    !*****************************************************************************************************
    
    subroutine calculate_H(H,p,q,u,ps,qs,us,gamma)
    
    parameter(NumEq = 8)
    real(8) H(NumEq),p(NumEq),q(NumEq),u(NumEq),gamma
    real(8) ps(NumEq),qs(NumEq),us(NumEq)
    real(8) pw(NumEq),qw(NumEq),uw(NumEq)
    
    uw = u - us
    pw = p - ps
    qw = q - qs
    
    !H(1) = u(1)*p(2) + u(2)*p(1) + u(1)*q(3) + u(3)*q(1)
    !H(2) = u(2)*p(2) + u(3)*q(2) + (u(7)*(p(7) - q(6)) + p(5))/u(1)
    !H(3) = u(2)*p(3) + u(3)*q(3) + (u(6)*(q(6) - p(7)) + q(5))/u(1)
    !H(4) = 0
    !H(5) = u(2)*p(5) + u(3)*q(5) + gamma*u(5)*(p(2) + q(3))
    !H(6) = u(3)*q(6) - u(2)*q(7) + u(6)*q(3) - u(7)*q(2)
    !H(7) = u(2)*p(7) - u(3)*p(6) + u(7)*p(2) - u(6)*p(3)
    !H(8) = 0
    
    H(1) = u(1)*pw(2) + uw(2)*p(1) + u(1)*qw(3) + uw(3)*q(1) + uw(1)*ps(2) + us(2)*pw(1) + uw(1)*qs(3) + us(3)*qw(1)
    H(2) = u(2)*pw(2) + u(3)*qw(2) + u(7)/u(1)*(pw(7) - qw(6)) + uw(2)*ps(2) + uw(3)*qs(2) + uw(7)/u(1)*(ps(7) - qs(6)) + uw(1)/u(1)*(us(2)*ps(2) + us(3)*qs(2)) + pw(5)/u(1)
    H(3) = u(2)*pw(3) + u(3)*qw(3) + u(6)/u(1)*(qw(6) - pw(7)) + uw(2)*ps(3) + uw(3)*qs(3) + uw(6)/u(1)*(qs(6) - ps(7)) + uw(1)/u(1)*(us(2)*ps(3) + us(3)*qs(3)) + qw(5)/u(1)
    H(4) = 0
    H(5) = u(2)*pw(5) + u(3)*qw(5) + gamma*u(5)*(pw(2) + qw(3)) + uw(2)*ps(5) + uw(3)*qs(5) + gamma*uw(5)*(ps(2) + qs(3))
    H(6) = u(3)*qw(6) - u(2)*qw(7) + uw(6)*q(3) - uw(7)*q(2) + uw(3)*qs(6) - uw(2)*qs(7) + us(6)*qw(3) - us(7)*qw(2)
    H(7) = u(2)*pw(7) - u(3)*pw(6) + uw(7)*p(2) - uw(6)*p(3) + uw(2)*ps(7) - uw(3)*ps(6) + us(7)*pw(2) - us(6)*pw(3)
    H(8) = 0
    
    end subroutine calculate_H
    
    !*****************************************************************************************************
    
    subroutine LF_Hamiltonian_2D
    
    use com
    
    Hhat1 = Hij - 0.5*max(abs(SR),abs(SL))*(p1wij - p2wij) - 0.5*max(abs(SU),abs(SD))*(q1wij - q2wij)
    
    end subroutine LF_Hamiltonian_2D
    
    !*****************************************************************************************************
    
    subroutine HLL_Hamiltonian_2D
    
    use com
    
    if (SL > 0 .and. SD > 0) then
        Hhat1 = HLD
    else if (SR < 0 .and. SU < 0) then
        Hhat1 = HRU
    else if (SL > 0 .and. SU < 0) then
        Hhat1 = HLU
    else if (SR < 0 .and. SD > 0) then
        Hhat1 = HRD
    else
        Hhat1 = (SL*SD*HRU - SL*SU*HRD - SR*SD*HLU + SR*SU*HLD)/((SR - SL)*(SU - SD)) + (SR*SL)/(SR - SL)*(p1wij - p2wij) + (SU*SD)/(SU - SD)*(q1wij - q2wij)
    end if
        
    end subroutine HLL_Hamiltonian_2D