    module com
    
    integer Nx,frameMAX,Lx,method
    real(8) pi
    real(8) xa,xb,t,dt,tend,CFL,maxu
    
    ! method = 1: DFT
    ! method = 2: FFT
    
    parameter(Lx = 16, frameMAX = 200, method = 1)
    parameter(Nx = 2*Lx - 1)
    
    parameter(pi = 4*atan(1.0d0))
    parameter(Nx1 = Nx + 1)

    real(8) X(0:Nx),hx,du(0:Nx)
    real(8) uh(0:Nx),dusin(Lx),ducos(Lx),uc
    real(8) sink(0:Nx,Lx),cosk(0:Nx,Lx)
    real(8) usin(Lx),ucos(Lx)
    real(8) Fu(0:Nx),Fusin(Lx),Fucos(Lx)
    real(8) ur(0:Nx),ui(0:Nx)
    
    end module com
    
    !*******************************
    
    module init1
    
    contains
    
    function u0(x)
    
    real(8) u0,x
    !u0 = 1 + cos(x) + cos(2*x) + cos(3*x) !+ cos(4*x)
    !u0 = 1 + sin(x) + sin(2*x) + sin(3*x)
    u0 = sin(cos(x))
    
    end function u0
    
    subroutine mesh
    
    use com
    
    xa = 0
    xb = 2*pi

    tend = 40
    
    end subroutine mesh
    
    end module init1
    
    !*******************************
    
    program main
    
    use com
    
    integer t1,t2
    
    call init_data
    
    call system_clock(t1)
    
    if (method == 1) then
        call Euler_Forward
    else if (method == 2) then
        call Euler_Forward2
    end if
    
    call system_clock(t2)
    
    print *,"time use is",(t2 - t1)/10000d0
    
    end program main

    !*******************************
    
    subroutine init_data
    
    use com
    
    use init1
    
    call mesh
    
    hx = (xb - xa)/Nx1
    
    uh = 0
    
    open(unit = 1,file = 'X.txt')
    
    do i = 0,Nx
        X(i) = xa + i*hx
        write(1,*) X(i)
    end do
    
    do i = 0,Nx
        uh(i) = u0(X(i))
    end do
    
    close(1)
    
    ! get_basis
    do i = 0,Nx
        do d = 1,Lx
            sink(i,d) = sin(d*X(i))
            cosk(i,d) = cos(d*X(i))
        end do
    end do
    
    call fft_u_to_uf(X,uh,Lx,ur,ui)
    
    usin = 0
    ucos = 0
    
    uc = ur(0)/Nx1
    usin(1:Lx - 1) = ui(1:Lx - 1)/(-Lx)
    ucos(1:Lx - 1) = ur(1:Lx - 1)/Lx
    
    uc = 0
    usin = 0
    ucos = 0
    ! L2 pro
    do i = 0,Nx
        uc = uc + uh(i)/Nx
        do d = 1,Lx
            usin(d) = usin(d) + uh(i)*sink(i,d)/Lx
            ucos(d) = ucos(d) + uh(i)*cosk(i,d)/Lx
        end do
    end do
    
    !call fft_u_to_uf(Xc,uh,Lx,ur,ui)
    !
    !do i = 0,Nx
    !    if (abs(ur(i)) < 1e-12) then
    !        ur(i) = 0
    !    end if
    !    if (abs(ui(i)) < 1e-12) then
    !        ui(i) = 0
    !    end if
    !end do
    
    !print *,"ur = "
    !print *,ur
    !print *,""
    !print *,"ui = "
    !print *,ui
    
    !call fft_uf_to_u(Xc,uh,Lx,ur,ui)
    
    uh = uc
    do d = 1,Lx
        uh = uh + usin(d)*sink(:,d) + ucos(d)*cosk(:,d)
    end do
    
    end subroutine init_data


    !*******************************   


    subroutine Lh
    
    use com
    
    uh = uc
    do d = 1,Lx
        uh = uh + usin(d)*sink(:,d) + ucos(d)*cosk(:,d)
    end do
    
    Fusin = 0
    Fucos = 0
    do d = 1,Lx
        do i = 0,Nx
            Fusin(d) = Fusin(d) + uh(i)*sink(i,d)/Lx
            Fucos(d) = Fucos(d) + uh(i)*cosk(i,d)/Lx
        end do
    end do
    
    do d = 1,Lx
        dusin(d) = -d*Fucos(d)
        ducos(d) = d*Fusin(d)
    end do
    dusin = -dusin
    ducos = -ducos
    
    end subroutine Lh


    !*******************************
    
    subroutine Euler_Forward
    
    use com
    
    real(8) t1
    
    CFL = 0.002
    dt = CFL*hx
    t = 0
    t1 = tend/frameMAX
    i1 = 1
    
    open(unit = 1,file = 'uh.txt')
    open(unit = 2,file = 'T.txt')
    
    do i = 0,Nx
        write(1,*) uh(i)
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
        
        usin = usin + dt*dusin
        ucos = ucos + dt*ducos
        
        call calculate_maxu
        
        print *,t,maxu
        
        if (t >= i1*t1) then
            do i = 0,Nx
                write(1,*) uh(i)
            end do
            write(2,*) t
            print *,"save the solution at t = ",t,i1*t1
            i1 = i1 + 1
        end if
        
    end do
    
    end subroutine Euler_Forward
    
    !*******************************

    subroutine Euler_Forward2
    
    use com
    
    real(8) t1
    
    CFL = 0.002
    dt = CFL*hx
    t = 0
    t1 = tend/frameMAX
    i1 = 1
    
    open(unit = 1,file = 'uh.txt')
    open(unit = 2,file = 'T.txt')
    
    do i = 0,Nx
        write(1,*) uh(i)
    end do
    write(2,*) t
    
    do while (t < tend)
        
        if (t + dt > tend) then
            dt = tend - t
            t = tend
        else
            t = t + dt
        end if
        
        call fftddx(du,uh)
        du = -du
        
        uh = uh + dt*du
        
        call calculate_maxu2
        
        print *,t,maxu
        
        if (t >= i1*t1) then
            do i = 0,Nx
                write(1,*) uh(i)
            end do
            write(2,*) t
            print *,"save the solution at t = ",t,i1*t1
            i1 = i1 + 1
        end if
        
    end do
    
    end subroutine Euler_Forward2

    !******************************* 

    subroutine calculate_maxu2
    
    use com
    
    maxu = 0
    
    do i = 0,Nx
            
        if (abs(uh(i)) > maxu) then
            maxu = abs(uh(i))
        end if
        
    end do
    
    end subroutine calculate_maxu2
    
    !******************************* 
    
    subroutine calculate_maxu
    
    use com
    
    maxu = 0
    
    uh = uc
    do d = 1,Lx
        uh = uh + usin(d)*sink(:,d) + ucos(d)*cosk(:,d)
    end do
    
    do i = 0,Nx
            
        if (abs(uh(i)) > maxu) then
            maxu = abs(uh(i))
        end if
        
    end do
    
    end subroutine calculate_maxu
    
    !******************************* 
    
    subroutine fft_u_to_uf(z,u,nl,ur,ui)
    
    dimension z(0:nl*2-1),u(0:nl*2-1)
    dimension ur(0:nl*2-1),ui(0:nl*2-1)
    double precision fin(0:nl*2-1),uff1(0:nl*2-1),uff2(0:nl*2-1)
    double complex fout(0:nl*2-1)
    integer iplan
    !      iplan=0
    nl2=nl*2

    !    ------------build fin----------

    do i=0,nl2-1
        fin(i)=u(i)
    end do

    !     --------------imkl fftw3-----------

    call dfftw_plan_dft_r2c_1d(iplan,nl2,fin,fout,64)
    call dfftw_execute_dft_r2c(iplan,fin,fout)
    call dfftw_destroy_plan(iplan)

    !     ---------conjugate------------

    !c      do 100 l=1,nl-1
    !c          fout(nl2-l)=dconjg(fout(l))
    !c100   continue

    !c     --------------separate real and imaginary---------------
    !c      do 12 l=0,nl2-1
    do l=0,nl
        ur(l)=real(fout(l))
        ui(l)=dimag(fout(l))
    end do

    return
    
    end subroutine fft_u_to_uf
    
    !******************************* 

    subroutine fft_uf_to_u(z,u,nl,ur,ui)
    
    dimension z(0:nl*2-1),u(0:nl*2-1)
    dimension ur(0:nl*2-1),ui(0:nl*2-1)
    double complex fin(0:nl*2-1)
    double precision fout(0:nl*2-1),uu(0:nl*2-1),dtemp
    integer iplan
    !      iplan=0
    nl2=nl*2

    !     --------build fin------------

    do l=0,nl
        fin(l)=dcmplx(ur(l),ui(l))
    end do

    !     -------------fftw3-------------

    call dfftw_plan_dft_c2r_1d(iplan,nl2,fin,fout,64)
    call dfftw_execute_dft_c2r(iplan,fin,fout)
    call dfftw_destroy_plan(iplan)

    !     ---------------rebuild u-------------

    do i=0,nl2-1
        u(i)=fout(i)/nl2
    end do

    return
    end subroutine fft_uf_to_u
    
    !******************************* 
    
    subroutine fftddx(dudx,uin)
    
    ! calculate the value of du/dx at grid point
    use com
    
    real(8) dudx(0:Nx),uin(0:Nx)
    real(8) ureal(0:Nx),uimag(0:Nx)
    real(8) uxreal(0:Nx),uximag(0:Nx)
    
    call fft_u_to_uf(X,uin,Lx,ureal,uimag)
    
    uxreal = 0
    uximag = 0
    do d = 1,Lx - 1
        uxreal(d) = -d*uimag(d)
        uximag(d) = d*ureal(d)
    end do
    
    call fft_uf_to_u(X,dudx,Lx,uxreal,uximag)
    
    end subroutine fftddx
    