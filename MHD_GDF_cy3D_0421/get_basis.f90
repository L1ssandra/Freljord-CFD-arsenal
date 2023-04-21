    subroutine get_basis
    
    include 'com.txt'
    
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
            phixG(i,j,2) = 1d0/hr1
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
            phiyG(i,j,3) = 1d0/hz1
            mm(3) = 1d0/3d0
                
            phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
            phiGLL(i,j,4,1) = lambdaL(i)**2 - 1d0/3d0
            phiGLL(i,j,4,2) = lambda(i)**2 - 1d0/3d0
            phiGR(j,4) = 2d0/3d0
            phiGL(j,4) = 2d0/3d0
            phiGU(i,4) = lambda(i)**2 - 1d0/3d0
            phiGD(i,4) = lambda(i)**2 - 1d0/3d0
            phixG(i,j,4) = 2d0*lambda(i)/hr1
            phiyG(i,j,4) = 0
            mm(4) = 4d0/45d0
                    
            phiG(i,j,5) = lambda(i)*lambda(j)
            phiGLL(i,j,5,1) = lambdaL(i)*lambda(j)
            phiGLL(i,j,5,2) = lambda(i)*lambdaL(j)
            phiGR(j,5) = lambda(j)
            phiGL(j,5) = -lambda(j)
            phiGU(i,5) = lambda(i)
            phiGD(i,5) = -lambda(i)
            phixG(i,j,5) = lambda(j)/hr1
            phiyG(i,j,5) = lambda(i)/hz1
            mm(5) = 1d0/9d0
                    
            phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,1) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,2) = lambdaL(j)**2 - 1d0/3d0
            phiGR(j,6) = lambda(j)**2 - 1d0/3d0
            phiGL(j,6) = lambda(j)**2 - 1d0/3d0
            phiGU(i,6) = 2d0/3d0
            phiGD(i,6) = 2d0/3d0
            phixG(i,j,6) = 0
            phiyG(i,j,6) = 2d0*lambda(j)/hz1
            mm(6) = 4d0/45d0
                
            phiG(i,j,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGLL(i,j,7,1) = lambdaL(i)**3 - 3d0*lambdaL(i)/5d0
            phiGLL(i,j,7,2) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGR(j,7) = 2d0/5d0
            phiGL(j,7) = -2d0/5d0
            phiGU(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGD(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phixG(i,j,7) = (3*lambda(i)**2 - 3d0/5d0)/hr1
            phiyG(i,j,7) = 0
            mm(7) = 4d0/175d0
                
            phiG(i,j,8) = (lambda(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,1) = (lambdaL(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,2) = (lambda(i)**2 - 1d0/3d0)*(lambdaL(j))
            phiGR(j,8) = (2d0/3d0)*(lambda(j))
            phiGL(j,8) = (2d0/3d0)*(lambda(j))
            phiGU(i,8) = (lambda(i)**2 - 1d0/3d0)
            phiGD(i,8) = -(lambda(i)**2 - 1d0/3d0)
            phixG(i,j,8) = 2d0*lambda(i)*lambda(j)/hr1
            phiyG(i,j,8) = (lambda(i)**2 - 1d0/3d0)/hz1
            mm(8) = 4d0/135d0
                
            phiG(i,j,9) = (lambda(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,1) = (lambdaL(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,2) = (lambda(i))*(lambdaL(j)**2 - 1d0/3d0)
            phiGR(j,9) = (lambda(j)**2 - 1d0/3d0)
            phiGL(j,9) = -(lambda(j)**2 - 1d0/3d0)
            phiGU(i,9) = lambda(i)*(2d0/3d0)
            phiGD(i,9) = lambda(i)*(2d0/3d0)
            phixG(i,j,9) = (lambda(j)**2 - 1d0/3d0)/hr1
            phiyG(i,j,9) = 2d0*lambda(i)*lambda(j)/hz1
            mm(9) = 4d0/135d0
                
            phiG(i,j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,1) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,2) = lambdaL(j)**3 - 3d0*lambdaL(j)/5d0
            phiGR(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGL(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGU(i,10) = 2d0/5d0
            phiGD(i,10) = -2d0/5d0
            phixG(i,j,10) = 0
            phiyG(i,j,10) = (3*lambda(j)**2 - 3d0/5d0)/hz1
            mm(10) = 4d0/175d0
        end do
        
        EzG(i,1) = 1
        EzG(i,2) = lambda(i)
        EzG(i,3) = lambda(i)**2 - 1d0/3d0
        
        EzxG(i,1) = 0
        EzxG(i,2) = 1/hr1
        EzxG(i,3) = 2*lambda(i)/hr1
        
        EzyG(i,1) = 0
        EzyG(i,2) = 1/hz1
        EzyG(i,3) = 2*lambda(i)/hz1
        
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
    
    do kk = 1,Lphi
        do i = 0,Nphi
            sink(i,kk) = sin(kk*Phi(i))
            cosk(i,kk) = cos(kk*Phi(i))
        end do
    end do
    
    do i1 = 1,NumGLP
        do j1 = 1,NumGLP
            do kk = 0,Nphi
                do n = 1,NumEq
                    RG(i1,j1,kk,n) = lambda(i1)
                end do
            end do
        end do
    end do
    
    do i = 0,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do j1 = 1,NumGLP
                    Redge1(i,j,kk,j1) = ra + i*hr
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                do i1 = 1,NumGLP
                    Redge2(i,j,kk,i1) = ra + (i - 0.5)*hr + hr1*lambda(i1)
                end do
            end do
        end do
    end do
    
    end subroutine get_basis