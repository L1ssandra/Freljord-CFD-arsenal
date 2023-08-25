    module com
    
    include 'mpif.h'
    
    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP,dimPk
    parameter(N_process = 16)
    parameter(Nx0 = 160, Ny0 = 160, Lphi = 0, kk = 3, NumEq = 2, NumGLP = 5, RKorder = 4, flux_type = 1)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)
    
    real(8) pi,gamma,gamma1
    parameter(dimPk = (kk + 1)*(kk + 2)/2)
    parameter(dimPk1 = (kk + 1)*(kk + 2)/2)
    parameter(Nphi = max(2*Lphi - 1,0))
    parameter(Nphi1 = Nphi + 1)
    !parameter(gamma = 1.4d0) ! other
    parameter(gamma = 5d0/3d0) ! jet, R-T
    parameter(gamma1 = gamma - 1)
    parameter(pi = 4*atan(1d0))
    
    ! The numerical solution and mesh
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,tRK,t1,t2,alphax,alphay,totaldiv,rij,eta
    real(8) uh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),du(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
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
    
    ! The Lh
    real(8) uGint3D(NumGLP,NumGLP,0:Nphi,NumEq),uGint(NumGLP,NumGLP,NumEq)
    real(8) pGint3D(NumGLP,NumGLP,0:Nphi,NumEq),qGint3D(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) mu1(NumGLP,NumGLP,0:Nphi,NumEq),mu2(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHSC(NumGLP,NumGLP,0:Nphi,NumEq),RHSCopen,RG(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHS(NumGLP,NumGLP,0:Nphi,NumEq),Fzsin(Lphi),Fzcos(Lphi),Fzzsin(Lphi),Fzzcos(Lphi)
    real(8) FR1(NumEq),FL1(NumEq),UR1(NumEq),UL1(NumEq),Fhat1(NumEq),SR,SL
    real(8) URstar(NumEq),ULstar(NumEq),Ustar(NumEq),UUstar(NumEq),UDstar(NumEq)
    real(8) URU1(NumEq),ULU1(NumEq),URD1(NumEq),ULD1(NumEq)
    real(8) URstarstar(NumEq),ULstarstar(NumEq),Ezhat
    real(8) EzVertex(0:Nx,0:Ny,0:Nphi)
    real(8) URU(0:Nx1,0:Ny1,0:Nphi,NumEq),ULU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URD(0:Nx1,0:Ny1,0:Nphi,NumEq),ULD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) L2(NumEq),L2pre(NumEq)
    
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
    
    ! sin
    module init1
    
    use com
    
    contains
    
    function B1(x,y,z)
    real(8) x,y,z,B1
    B1 = -sin(y)
    end function B1
    
    function B2(x,y,z)
    real(8) x,y,z,B2
    B2 = sin(x)
    end function B2
    
    subroutine mesh
    
    use com
    
    eta = 0.05
    
    xa = 0
    xb = 2*pi
    ya = 0
    yb = 2*pi

    bcR = 1
    bcL = 1
    bcU = 1
    bcD = 1

    tend = 5
    
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
    
    ! The initial value:
    ! 1: sin
    
    use init1
    
    real(8) U1
    U1(x,y,z) = B1(x,y,z)
    real(8) U2
    U2(x,y,z) = B2(x,y,z)
    
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
    
    ! L2 Pro for Uh
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh(i,j,k,d,1) = uh(i,j,k,d,1) + 0.25*weight(i1)*weight(j1)*U1(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,2) = uh(i,j,k,d,2) + 0.25*weight(i1)*weight(j1)*U2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        uh(:,:,:,d,:) = uh(:,:,:,d,:)/mm(d)
    end do
    
    end subroutine init_data
    
    !*****************************************************************************************************
    
    subroutine save_solution
    
    use com
    
    real uhsave(NumEq)
    integer the_idx1,the_idy1
    
    !uh = du
    
    if (myid1 == 1) then
        
        open(unit = 1,file = 'Bx.txt')
        open(unit = 2,file = 'By.txt')
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
                end if
                    
            end do
        end do
    !end do
    
    if (myid1 == 1) then
        close(1)
        close(2)
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
        end do
        
    end do
    
    end subroutine get_basis
    
    !*****************************************************************************************************
    
    subroutine calculate_L2_Error
    
    use com
    
    use init1
    
    real(8) U1
    U1(x,y,z) = B1(x,y,z)
    real(8) U2
    U2(x,y,z) = B2(x,y,z)
    
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
                        L2(1) = L2(1) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,1) - exp(-eta*tend)*U1(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                        L2(2) = L2(2) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,2) - exp(-eta*tend)*U2(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
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
        print *,"Bx : ",L2(1)
        print *,"By : ",L2(2)
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
    
    CFL = 0.2
    t = 0
    
    call div_free
    
    call calculate_umax
        
    if (myid1 == 1) then
        open(unit = 12,file = 'Latest_result.txt')
        print *,t,umax
        write(12,*) t,umax
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
        uII = uh
        
        ! Stage I
        do i = 1,5
            
            call Lh
            
            uI = uh + (dt/6d0)*du
            
            tRK = tRK + (dt/6d0)
            
            uh = uI
            
            call div_free
            
        end do
        
        uII = 0.04d0*uII + 0.36d0*uI
        
        uI = 15*uII - 5*uI
        
        uh = uI
        
        tRK = tRK - 0.5*dt
        
        ! Stage II
        do i = 6,9
            
            call Lh
            
            uI = uh + (dt/6d0)*du
            
            tRK = tRK + dt/6d0
            
            uh = uI
            
            call div_free
            
        end do
        
        call Lh
        
        uh = uII + 0.6d0*uI + (dt/10d0)*du
        
        call div_free
        
        call calculate_umax
        
        if (myid1 == 1) then
            print *,t,umax
            write(12,*) t,umax
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
    
    alphax = 1d0
    alphay = 1d0
    
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
    
    end subroutine set_bc
    
    !*****************************************************************************************************
    
    subroutine Lh
    
    use com
    
    real(8) Fx(NumGLP,NumGLP,0:Nphi,NumEq), Fy(NumGLP,NumGLP,0:Nphi,NumEq), Fz(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) rhoij,uij,vij,wij,Eij,B1ij,B2ij,B3ij,pij,Sij,Tij,Kij,rB1ij,rB2ij,rB3ij
    
    real(8),allocatable :: UR(:,:,:,:,:),UL(:,:,:,:,:),UU(:,:,:,:,:),UD(:,:,:,:,:)
    real(8),allocatable :: pR(:,:,:,:,:),pL(:,:,:,:,:),qU(:,:,:,:,:),qD(:,:,:,:,:)
    real(8),allocatable :: pU(:,:,:,:,:),pD(:,:,:,:,:),qR(:,:,:,:,:),qL(:,:,:,:,:)
    real(8),allocatable :: FR(:,:,:,:,:),FL(:,:,:,:,:),FU(:,:,:,:,:),FD(:,:,:,:,:)
    real(8),allocatable :: mu1R(:,:,:,:,:),mu1L(:,:,:,:,:),mu2U(:,:,:,:,:),mu2D(:,:,:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:),mu1hat(:,:,:,:,:),mu2hat(:,:,:,:,:)
    real(8),allocatable :: uhxhat(:,:,:,:,:),uhyhat(:,:,:,:,:),phxhat(:,:,:,:,:),qhyhat(:,:,:,:,:)
    real(8),allocatable :: ph(:,:,:,:,:),qh(:,:,:,:,:),uhpre(:,:,:,:,:)
    
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
    
    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(uhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(phxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(qhyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu1hat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(mu2hat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    
    call set_bc
    
    uij = 1
    vij = 1
    
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
                            qh(i,j,k,d,n) = qh(i,j,k,d,n) + (0.5d0/hy)*weight(i1)*(uhyhat(i,j,k,i1,n)*phiGU(i1,d) - uhyhat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        ph(:,:,:,d,:) = ph(:,:,:,d,:)/mm(d)
        qh(:,:,:,d,:) = qh(:,:,:,d,:)/mm(d)
    end do
    
    uhpre = uh
    
    uh = ph
    call set_bc
    ph = uh
    
    uh = qh
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
            do n = 1,NumEq
                do d = 1,dimPk
                    do k = 0,Nphi
                        uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                        pGint3D(:,:,k,n) = pGint3D(:,:,k,n) + ph(i,j,k,d,n)*phiG(:,:,d)
                        qGint3D(:,:,k,n) = qGint3D(:,:,k,n) + qh(i,j,k,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            do k = 0,Nphi
                do j1 = 1,NumGLP
                    do i1 = 1,NumGLP
                        B1ij = uGint3D(i1,j1,k,1)
                        B2ij = uGint3D(i1,j1,k,2)
                        
                        Fx(i1,j1,k,1) = 0
                        Fx(i1,j1,k,2) = -(vij*B1ij - uij*B2ij)
                        
                        Fy(i1,j1,k,1) = (vij*B1ij - uij*B2ij)
                        Fy(i1,j1,k,2) = 0
                        
                        mu1(i1,j1,k,1) = 0
                        mu1(i1,j1,k,2) = pGint3D(i1,j1,k,2)*eta - qGint3D(i1,j1,k,1)*eta
                        
                        mu2(i1,j1,k,1) = qGint3D(i1,j1,k,1)*eta - pGint3D(i1,j1,k,2)*eta
                        mu2(i1,j1,k,2) = 0
                    end do
                end do
            end do
            
            do n = 1,NumEq
                do d = 1,dimPk1
                    do k = 0,Nphi
                        do j1 = 1,NumGLP
                            do i1 = 1,NumGLP
                                if (d > 1) then
                                    du(i,j,k,d,n) = du(i,j,k,d,n) + 0.25d0*weight(i1)*weight(j1)*((Fx(i1,j1,k,n) - mu1(i1,j1,k,n))*phixG(i1,j1,d) + (Fy(i1,j1,k,n) - mu2(i1,j1,k,n))*phiyG(i1,j1,d))
                                end if
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
                    B1ij = uR(i,j,k,j1,1)
                    B2ij = uR(i,j,k,j1,2)
                    FR(i,j,k,j1,1) = 0
                    FR(i,j,k,j1,2) = -(vij*B1ij - uij*B2ij)
                    mu1R(i,j,k,j1,1) = 0
                    mu1R(i,j,k,j1,2) = pR(i,j,k,j1,2)*eta - qR(i,j,k,j1,1)*eta
                end do
            end do
        end do
    end do
    
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 1,Nx1
                    B1ij = uL(i,j,k,j1,1)
                    B2ij = uL(i,j,k,j1,2)
                    FL(i,j,k,j1,1) = 0
                    FL(i,j,k,j1,2) = -(vij*B1ij - uij*B2ij)
                    mu1L(i,j,k,j1,1) = 0
                    mu1L(i,j,k,j1,2) = pL(i,j,k,j1,2)*eta - qL(i,j,k,j1,1)*eta
                end do
            end do
        end do
    end do
    
    ! The y-Flux
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    B1ij = UU(i,j,k,i1,1)
                    B2ij = UU(i,j,k,i1,2)
                    FU(i,j,k,i1,1) = vij*B1ij - uij*B2ij
                    FU(i,j,k,i1,2) = 0
                    mu2U(i,j,k,i1,1) = qU(i,j,k,i1,1)*eta - pU(i,j,k,i1,2)*eta
                    mu2U(i,j,k,i1,2) = 0
                end do
            end do
        end do
    end do
    
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny1
                do i = 1,Nx
                    B1ij = UD(i,j,k,i1,1)
                    B2ij = UD(i,j,k,i1,2)
                    FD(i,j,k,i1,1) = vij*B1ij - uij*B2ij
                    FD(i,j,k,i1,2) = 0
                    mu2D(i,j,k,i1,1) = qD(i,j,k,i1,1)*eta - pD(i,j,k,i1,2)*eta
                    mu2D(i,j,k,i1,2) = 0
                end do
            end do
        end do
    end do
    
    ! calculate Fx hat
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    SR = 1
                    SL = -1
                    FR1 = FL(i + 1,j,k,j1,:)
                    FL1 = FR(i,j,k,j1,:)
                    UR1 = UL(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    end if
                    Fxhat(i,j,k,j1,:) = Fhat1
                    mu1hat(i,j,k,j1,:) = 0.5*(mu1R(i,j,k,j1,:) + mu1L(i + 1,j,k,j1,:))
                    !phxhat(i,j,k,j1,:) = 0.5*(pR(i,j,k,j1,:) + pL(i + 1,j,k,j1,:))
                end do
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    SR = 1
                    SL = -1
                    FR1 = FD(i,j + 1,k,i1,:)
                    FL1 = FU(i,j,k,i1,:)
                    UR1 = UD(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    end if
                    Fyhat(i,j,k,i1,:) = Fhat1
                    mu2hat(i,j,k,i1,:) = 0.5*(mu2U(i,j,k,i1,:) + mu2D(i,j + 1,k,i1,:))
                    !qhyhat(i,j,k,i1,:) = 0.5*(qU(i,j,k,i1,:) + qD(i,j + 1,k,i1,:))
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
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hx)*weight(j1)*((Fxhat(i,j,k,j1,n) - mu1hat(i,j,k,j1,n))*phiGR(j1,d) - (Fxhat(i - 1,j,k,j1,n) - mu1hat(i - 1,j,k,j1,n))*phiGL(j1,d))
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
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hy)*weight(i1)*((Fyhat(i,j,k,i1,n) - mu2hat(i,j,k,i1,n))*phiGU(i1,d) - (Fyhat(i,j - 1,k,i1,n) - mu2hat(i,j - 1,k,i1,n))*phiGD(i1,d))
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
    