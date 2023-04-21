    subroutine div_free_Balsara
    
    include 'com.txt'
    
    real BxR(k + 1),BxL(k + 1),ByU(k + 1),ByD(k + 1)
    real Bxint(dimPk),Byint(dimPk)
    real a0R,a1R,a2R,a0L,a1L,a2L
    real b0U,b1U,b2U,b0D,b1D,b2D
    real a00,a10,a01,a20,a11,a02,a30,a21,a12
    real b00,b10,b01,b20,b11,b02,b21,b12,b03
    real rxy,ryx,ri,S1,S2
    
    uh(:,:,:,:,6:7) = 0
    
    do i = 1,Nx
        do j = 1,Ny
            
            divphi = 0
            Bphiphi = 0
            uGint3D = 0
            do kk = 0,Nphi
                do d = 1,dimPk
                    uGint3D(:,:,kk,8) = uGint3D(:,:,kk,8) + uh(i,j,kk,d,8)*phiG(:,:,d)
                end do
            end do
            
            ! Fourier spectral method to calculate d(Bphi)/d(phi)
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                        
                    Fzsin = 0
                    Fzcos = 0
                    Fzzsin = 0
                    Fzzcos = 0
                    do d = 1,Lphi
                        do kk = 0,Nphi
                            Fzsin(d) = Fzsin(d) + sink(kk,d)*uGint3D(i1,j1,kk,8)
                            Fzcos(d) = Fzcos(d) + cosk(kk,d)*uGint3D(i1,j1,kk,8)
                        end do
                    end do
                    Fzsin = Fzsin/Lphi
                    Fzcos = Fzcos/Lphi
                    
                    do d = 1,Lphi
                        Fzzsin(d) = -d*Fzcos(d)
                        Fzzcos(d) = d*Fzsin(d)
                    end do
                        
                    do kk = 0,Nphi
                        do d = 1,Lphi
                            divphi(i1,j1,kk) = divphi(i1,j1,kk) + Fzzsin(d)*sink(kk,d) + Fzzcos(d)*cosk(kk,d)
                        end do
                    end do
                    
                end do
            end do
            
            ! project d(Bphi)/d(phi) to (r,z) plane
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            Bphiphi(kk,d) = Bphiphi(kk,d) + 0.25*weight(i1)*weight(j1)*divphi(i1,j1,kk)*phiG(i1,j1,d)/mm(d)
                        end do
                    end do
                end do
            end do
            
            do kk = 0,Nphi
            
                BxR = Br(i,j,kk,:)
                BxL = Br(i - 1,j,kk,:)
                ByU = Bz(i,j,kk,:)
                ByD = Bz(i,j - 1,kk,:)
            
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
                
                c00 = -Bphiphi(kk,1)
                c10 = -Bphiphi(kk,2)
                c01 = -Bphiphi(kk,3)
                c20 = -Bphiphi(kk,4)
                c11 = -Bphiphi(kk,5)
                c02 = -Bphiphi(kk,6)
            
                rxy = hr/hz
                ryx = hz/hr
                ri = Rc(i)
    
                ! The reconstruction of B = (Bx,By) from the interface
                a11 = (1d0/2d0)*(a1R - a1L)
                a02 = (1d0/2d0)*(a2R + a2L)
                a12 = (1d0/2d0)*(a2R - a2L)
                
                b11 = (1d0/2d0)*(b1U - b1D)
                b20 = (1d0/2d0)*(b2U + b2D)
                b21 = (1d0/2d0)*(b2U - b2D)
                
                a20 = (hr/4)*(c10 - 2d0/hz*b11)
                a21 = (hr/8)*c11
                b12 = (hz/8)*c11
                a30 = (hr/6)*(c20 - 2d0/hz*b21)
                b02 = (hz/4)*(c01 - 2d0/hr*a11)
                b03 = (hz/6)*(c02 - 2d0/hr*a12)
                a01 = (1d0/2d0)*(a1R + a1L) - (2d0/3d0)*a21
                b10 = (1d0/2d0)*(b1U + b1D) - (2d0/3d0)*b12
                b01 = (1d0/2d0)*(b0U - b0D) - (2d0/5d0)*b03
                b00 = (1d0/2d0)*(b0U + b0D) - (2d0/3d0)*b02
                a10 = (1d0/2d0)*(a0R - a0L) - (2d0/5d0)*a30
                a00 = (1d0/2d0)*(a0R + a0L) - (2d0/3d0)*a20
                
                uh(i,j,kk,1,6) = a00
                uh(i,j,kk,2,6) = a10
                uh(i,j,kk,3,6) = a01
                uh(i,j,kk,4,6) = a20
                uh(i,j,kk,5,6) = a11
                uh(i,j,kk,6,6) = a02
                uh(i,j,kk,7,6) = a30
                uh(i,j,kk,8,6) = a21
                uh(i,j,kk,9,6) = a12
            
                uh(i,j,kk,1,7) = b00
                uh(i,j,kk,2,7) = b10
                uh(i,j,kk,3,7) = b01
                uh(i,j,kk,4,7) = b20
                uh(i,j,kk,5,7) = b11
                uh(i,j,kk,6,7) = b02
                uh(i,j,kk,8,7) = b21
                uh(i,j,kk,9,7) = b12
                uh(i,j,kk,10,7) = b03
                
            end do
        end do
    end do
    
    
    end subroutine div_free_Balsara