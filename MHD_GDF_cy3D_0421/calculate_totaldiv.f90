    subroutine calculate_totaldiv
    
    include 'com.txt'
    
    totaldiv = 0
    open(unit = 10,file = 'div.txt')
    
    do j = 1,Nz
        do i = 1,Nr
            
            uGdiv = 0
            divphi = 0
            
            do kk = 0,Nphi
                do d = 1,dimPk
                    uGdiv(:,:,kk) = uGdiv(:,:,kk) + uh(i,j,kk,d,6)*phixG(:,:,d) + uh(i,j,kk,d,7)*phiyG(:,:,d)
                end do
            end do
            
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
            
            uGdiv = (uGdiv + divphi)/(Rc(i) + hr1*RG(:,:,:,1))
            
            div0 = 0
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do kk = 0,Nphi
                        totaldiv = totaldiv + weight(i1)*weight(j1)*uGdiv(i1,j1,kk)**2
                        div0 = div0 + weight(i1)*weight(j1)*uGdiv(i1,j1,kk)**2
                    end do
                end do
            end do
            
            write(10,*) sqrt(0.25*div0)
            
        end do
    end do
    
    close(10)
    totaldiv = (0.25*totaldiv/(Nx*Ny*Nphi1))**0.5
    
    end subroutine calculate_totaldiv