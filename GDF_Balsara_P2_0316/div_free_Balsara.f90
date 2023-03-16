    subroutine div_free_Balsara
    
    include 'com.txt'
    
    real BxR(k + 1),BxL(k + 1),ByU(k + 1),ByD(k + 1)
    real Bxint(dimPk),Byint(dimPk)
    real a0R,a1R,a2R,a0L,a1L,a2L
    real b0U,b1U,b2U,b0D,b1D,b2D
    real a00,a10,a01,a20,a11,a02,a30,a12
    real b00,b10,b01,b20,b11,b02,b21,b03
    real rxy,ryx
    
    uh(:,:,:,6:7) = 0
    
    do i = 1,Nx
        do j = 1,Ny
            
            BxR = Bx(i,j,:)
            BxL = Bx(i - 1,j,:)
            ByU = By(i,j,:)
            ByD = By(i,j - 1,:)
            
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
            
            uh(i,j,1,6) = a00
            uh(i,j,2,6) = a10
            uh(i,j,3,6) = a01
            uh(i,j,4,6) = a20
            uh(i,j,5,6) = a11
            uh(i,j,6,6) = a02
            uh(i,j,7,6) = a30
            uh(i,j,9,6) = a12
            
            uh(i,j,1,7) = b00
            uh(i,j,2,7) = b10
            uh(i,j,3,7) = b01
            uh(i,j,4,7) = b20
            uh(i,j,5,7) = b11
            uh(i,j,6,7) = b02
            uh(i,j,8,7) = b21
            uh(i,j,10,7) = b03
        end do
    end do
    
    
    end subroutine div_free_Balsara