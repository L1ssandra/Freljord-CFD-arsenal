    subroutine LhEz
    
    include 'com.txt'
    
    ! calculate the flux Ez(u-,u+)
    do i = 0,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do j1 = 1,NumGLP
                    SRmax = UR(i,j,kk,j1,2)/UR(i,j,kk,j1,1)
                    SLmax = UL(i + 1,j,kk,j1,2)/UL(i + 1,j,kk,j1,1)
                    ! L-F Flux
                    SR = max(abs(SRmax),abs(SLmax))
                    EphiRL(i,j,kk,j1) = -0.5*(FR(i,j,kk,j1,7) + FL(i + 1,j,kk,j1,7) - SR*(UL(i + 1,j,kk,j1,7) - UR(i,j,kk,j1,7)))
                    EzRL(i,j,kk,j1) = 0.5*(FR(i,j,kk,j1,8) + FL(i + 1,j,kk,j1,8) - SR*(UL(i + 1,j,kk,j1,8) - UR(i,j,kk,j1,8)))
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                do i1 = 1,NumGLP
                    SRmax = UU(i,j,kk,i1,3)/UU(i,j,kk,i1,1)
                    SLmax = UD(i,j + 1,kk,i1,3)/UD(i,j + 1,kk,i1,1)
                    ! L-F Flux
                    SR = max(abs(SRmax),abs(SLmax))
                    EphiUD(i,j,kk,i1) = 0.5*(FU(i,j,kk,i1,6) + FD(i,j + 1,kk,i1,6) - SR*(UD(i,j + 1,kk,i1,6) - UU(i,j,kk,i1,6)))
                    ErUD(i,j,kk,i1) = -0.5*(FU(i,j,kk,i1,8) + FD(i,j + 1,kk,i1,8) - SR*(UD(i,j + 1,kk,i1,8) - UU(i,j,kk,i1,8)))
                end do
            end do
        end do
    end do
    
    ! calculate d(Bphi)
    du(:,:,:,:,8) = du8
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do j1 = 1,NumGLP
                        du(i,j,kk,d,8) = du(i,j,kk,d,8) - (0.5d0/hr)*weight(j1)*(EzRL(i,j,kk,j1)*phiGR(j1,d) - EzRL(i - 1,j,kk,j1)*phiGL(j1,d))
                    end do
                    do i1 = 1,NumGLP
                        du(i,j,kk,d,8) = du(i,j,kk,d,8) + (0.5d0/hz)*weight(i1)*(ErUD(i,j,kk,i1)*phiGU(i1,d) - ErUD(i,j - 1,kk,i1)*phiGD(i1,d))
                    end do
                end do
            end do
        end do
    end do
    
    do d = 1,dimPk1
        du(:,:,:,d,8) = du(:,:,:,d,8)/mm(d)
    end do
    
    !EphiRL = -Fxhat(:,:,:,:,7)
    !EphiUD = Fyhat(:,:,:,:,6)
    
    !EzRL = Fxhat(:,:,:,:,8)
    !ErUD = -Fyhat(:,:,:,:,8)
    
    RHS1 = 0
    RHS2 = 0
    
    ! calculate d(E)/d(phi)
    do i = 0,Nx
        do j = 1,Ny
            do j1 = 1,NumGLP
                Fzsin = 0
                Fzcos = 0
                Fzzsin = 0
                Fzzcos = 0
                do d = 1,Lphi
                    do kk = 0,Nphi
                        Fzsin(d) = Fzsin(d) + sink(kk,d)*EzRL(i,j,kk,j1)
                        Fzcos(d) = Fzcos(d) + cosk(kk,d)*EzRL(i,j,kk,j1)
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
                        RHS1(i,j,kk,j1) = RHS1(i,j,kk,j1) + Fzzsin(d)*sink(kk,d) + Fzzcos(d)*cosk(kk,d)
                    end do
                end do
                        
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do i1 = 1,NumGLP
                Fzsin = 0
                Fzcos = 0
                Fzzsin = 0
                Fzzcos = 0
                do d = 1,Lphi
                    do kk = 0,Nphi
                        Fzsin(d) = Fzsin(d) + sink(kk,d)*ErUD(i,j,kk,i1)
                        Fzcos(d) = Fzcos(d) + cosk(kk,d)*ErUD(i,j,kk,i1)
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
                        RHS2(i,j,kk,i1) = RHS2(i,j,kk,i1) - Fzzsin(d)*sink(kk,d) - Fzzcos(d)*cosk(kk,d)
                    end do
                end do
                        
            end do
        end do
    end do
    
    ! calculate each Uh at vertex
    URU = 0
    ULU = 0
    URD = 0
    ULD = 0
    do i = 0,Nx1
        do j = 0,Ny1
            do kk = 0,Nphi
                do d = 1,dimPk
                    URU(i,j,kk,:) = URU(i,j,kk,:) + uh(i,j,kk,d,:)*phiRU(d)
                    ULU(i,j,kk,:) = ULU(i,j,kk,:) + uh(i,j,kk,d,:)*phiLU(d)
                    URD(i,j,kk,:) = URD(i,j,kk,:) + uh(i,j,kk,d,:)*phiRD(d)
                    ULD(i,j,kk,:) = ULD(i,j,kk,:) + uh(i,j,kk,d,:)*phiLD(d)
                end do
            end do
        end do
    end do
    
    ! calculate Ez at vertex
    do i = 0,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                URU1 = ULD(i + 1,j + 1,kk,:)
                ULU1 = URD(i,j + 1,kk,:)
                URD1 = ULU(i + 1,j,kk,:)
                ULD1 = URU(i,j,kk,:)
            
                call LF_Flux_2D
            
                EzVertex(i,j,kk) = Ezhat
            end do
        end do
    end do
    
    dBr = 0
    dBz = 0
    
    ! DG scheme of Ez1
    do i = 0,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do d = 1,k + 1
                    do j1 = 1,NumGLP
                        dBr(i,j,kk,d) = dBr(i,j,kk,d) + 0.5*weight(j1)*EphiRL(i,j,kk,j1)*EzyG(j1,d)
                        dBr(i,j,kk,d) = dBr(i,j,kk,d) + 0.5*weight(j1)*RHS1(i,j,kk,j1)*EzG(j1,d)
                        if (RHSCopen == 1) then
                            dBr(i,j,kk,d) = dBr(i,j,kk,d) + 0.5*weight(j1)*(ra + i*hr)*S6(ra + i*hr,Zc(j) + hz1*lambda(j1),Phi(kk),tRK)*EzG(j1,d)
                        end if
                    end do
                    dBr(i,j,kk,d) = (dBr(i,j,kk,d) - EzVertex(i,j,kk)*EzU(d)/hz + EzVertex(i,j - 1,kk)*EzD(d)/hz)/mmE(d)
                end do
            end do
        end do
    end do
    
    ! DG scheme of Ez2
    do i = 1,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                do d = 1,k + 1
                    do i1 = 1,NumGLP
                        dBz(i,j,kk,d) = dBz(i,j,kk,d) - 0.5*weight(i1)*EphiUD(i,j,kk,i1)*EzxG(i1,d)
                        dBz(i,j,kk,d) = dBz(i,j,kk,d) + 0.5*weight(i1)*RHS2(i,j,kk,i1)*EzG(i1,d)
                        if (RHSCopen == 1) then
                            dBz(i,j,kk,d) = dBz(i,j,kk,d) + 0.5*weight(i1)*(Rc(i) + hr1*lambda(i1))*S7(Rc(i) + hr1*lambda(i1),za + j*hz,Phi(kk),tRK)*EzG(i1,d)
                        end if
                    end do
                    dBz(i,j,kk,d) = (dBz(i,j,kk,d) + EzVertex(i,j,kk)*EzR(d)/hr - EzVertex(i - 1,j,kk)*EzL(d)/hr)/mmE(d)
                end do
            end do
        end do
    end do
    
    end subroutine LhEz
    
    