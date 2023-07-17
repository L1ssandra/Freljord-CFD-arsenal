    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP
    parameter(N_process = 16)
    parameter(Nx0 = 512, Ny0 = 512, Lphi = 0, kk = 2, NumEq = 8, NumGLP = 3, RKorder = 4, flux_type = 2)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) pi
    parameter(dimPk = (kk + 2)*(kk + 3)/2)
    parameter(dimPk1 = (kk + 2)*(kk + 3)/2)
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
    
    use init1
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
    
    ! L2 Pro
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
    
    uh(1:Nx,1:Ny,:,:,:) = uh0((myidx - 1)*Nx + 1:myidx*Nx,(myidy - 1)*Ny + 1:myidy*Ny,:,:,:)
    
    end subroutine init_data
    
    !*****************************************************************************************************
    
    subroutine save_solution
    
    use com
    
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
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
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
    
    CFL = 0.95
    t = 0
    
    do while (t < tend)
        
        call calculate_dt
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        uI = uh
        uII = uh
        
        ! Stage I
        do i = 1,5
            uh = uI
            call Lh
            uI = uh + (dt/6d0)*du
        end do
        
        uII = 0.04d0*uII + 0.36d0*uI
        uI = 15*uII - 5*uI
        
        ! Stage II
        do i = 6,9
            uh = uI
            call Lh
            uI = uh + (dt/6d0)*du
        end do
        
        uh = uI
        call Lh
        uh = uII + 0.6d0*uI + (dt/10d0)*du
        
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
                                    call MPI_SEND(uh(1,1:Ny,k,d,n),Ny,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                                end if
                            end if
            
                        else
                
                            the_idx = i + 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1,1:Ny,k,d,n),Ny,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(Nx1,1:Ny,k,d,n),Ny,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                        end if
                        
                        ! The Left condition
                        if (i == 1) then
                            the_idx = Nx_process
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcL == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(Nx,1:Ny,k,d,n),Ny,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i - 1
                            the_idy = j
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(Nx,1:Ny,k,d,n),Ny,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0,1:Ny,k,d,n),Ny,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                        end if
            
                        ! The Up condition
                        if (j == Ny_process) then
                            the_idx = i
                            the_idy = 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcU == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(1:Nx,1,k,d,n),Nx,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j + 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1:Nx,1,k,d,n),Nx,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(1:Nx,Ny1,k,d,n),Nx,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                        end if
            
                        ! The Down condition
                        if (j == 1) then
                            the_idx = i
                            the_idy = Ny_process
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (bcD == 1) then
                                if (myid1 == the_id2) then
                                    call MPI_SEND(uh(1:Nx,Ny,k,d,n),Nx,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                                end if
                            end if
                        else
                            the_idx = i
                            the_idy = j - 1
                            the_id2 = the_idx + Nx_process*(the_idy - 1)
                
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1:Nx,Ny,k,d,n),Nx,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                            end if
                
                        end if
            
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(1:Nx,0,k,d,n),Nx,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
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
    real(8),allocatable :: URU(:,:,:,:),ULU(:,:,:,:),URD(:,:,:,:),ULD(:,:,:,:)
    real(8),allocatable :: EzVertex(:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:)
    real(8),allocatable :: EzR1(:,:,:,:),EzL1(:,:,:,:),ErU(:,:,:,:),ErD(:,:,:,:)
    real(8),allocatable :: EphiRL(:,:,:,:), EzRL(:,:,:,:), EphiUD(:,:,:,:), ErUD(:,:,:,:)
    real(8),allocatable :: RHS1(:,:,:,:), RHS2(:,:,:,:)
    
    allocate(UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(UD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(FR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(FD(Nx,Ny1,0:Nphi,NumGLP,NumEq))
    
    allocate(URU(0:Nx1,0:Ny1, 0:Nphi, NumEq))
    allocate(ULU(0:Nx1,0:Ny1,0:Nphi,NumEq))
    allocate(URD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    allocate(ULD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    allocate(EzVertex(0:Nx,0:Ny,0:Nphi))
    
    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    
    allocate(EzR1(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(EzL1(1:Nx1,Ny,0:Nphi,NumGLP))
    allocate(ErU(Nx,0:Ny,0:Nphi,NumGLP))
    allocate(ErD(Nx,0:Ny,0:Nphi,NumGLP))
    
    allocate(EphiRL(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(EzRL(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(EphiUD(Nx,0:Ny,0:Nphi,NumGLP))
    allocate(ErUD(Nx,0:Ny,0:Nphi,NumGLP))
    
    allocate(RHS1(0:Nx,Ny,0:Nphi,NumGLP))
    allocate(RHS2(Nx,0:Ny,0:Nphi,NumGLP))
    
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
                        !call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 1
                        !call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 1
                        !call HLLD_Flux
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
                        !call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 2
                        !call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 2
                        !call HLLD_Flux
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