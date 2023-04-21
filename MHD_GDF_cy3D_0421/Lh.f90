    subroutine Lh
    
    include 'com.txt'
    
    du = 0
    
    !call real_solution
    
    ! The value of num solution on the GL points
    
    !$omp parallel default(shared) private(uGint3D,RHS,RHSC,rhoM,uM,vM,wM,EM,B1M,B2M,B3M,pM,SM,TM,KM,Fx,Fy,Fz,i,j,d,n,i1,j1,kk,Fzsin,Fzcos,Fzzsin,Fzzcos,rB1M,rB2M,rB3M)
    
    !$omp do collapse(2)
    do i = 1,Nx
        do j = 1,Ny
            
            uGint3D = 0
            RHS = 0
            do kk = 0,Nphi
                do n = 1,NumEq
                    do d = 1,dimPk
                        uGint3D(:,:,kk,n) = uGint3D(:,:,kk,n) + uh(i,j,kk,d,n)*phiG(:,:,d)
                    end do
                end do
            end do
            
            rhoM = uGint3D(:,:,:,1)
            uM = uGint3D(:,:,:,2)/rhoM
            vM = uGint3D(:,:,:,3)/rhoM
            wM = uGint3D(:,:,:,4)/rhoM
            EM = uGint3D(:,:,:,5)
            rB1M = uGint3D(:,:,:,6)
            rB2M = uGint3D(:,:,:,7)
            B1M = rB1M/(Rc(i) + hr1*RG(:,:,:,1))
            B2M = rB2M/(Rc(i) + hr1*RG(:,:,:,1))
            B3M = uGint3D(:,:,:,8)
            rB3M = B3M*(Rc(i) + hr1*RG(:,:,:,1))
    
            pM = gamma1*(EM - 0.5d0*rhoM*(uM**2 + vM**2 + wM**2) - 0.5d0*(B1M**2 + B2M**2 + B3M**2))
    
            SM = pM + 0.5d0*(B1M**2 + B2M**2 + B3M**2)
            TM = EM + SM
            KM = uM*B1M + vM*B2M + wM*B3M
    
            Fx(:,:,:,1) = uGint3D(:,:,:,2)
            Fx(:,:,:,2) = rhoM*uM**2 + SM - B1M**2
            Fx(:,:,:,3) = rhoM*uM*vM - B1M*B2M
            Fx(:,:,:,4) = rhoM*uM*wM - B1M*B3M
            Fx(:,:,:,5) = TM*uM - KM*B1M
            Fx(:,:,:,6) = 0
            Fx(:,:,:,7) = uM*rB2M - vM*rB1M
            Fx(:,:,:,8) = uM*B3M - wM*B1M
    
            Fy(:,:,:,1) = uGint3D(:,:,:,3)
            Fy(:,:,:,2) = rhoM*uM*vM - B1M*B2M
            Fy(:,:,:,3) = rhoM*vM**2 + SM - B2M**2
            Fy(:,:,:,4) = rhoM*vM*wM - B2M*B3M
            Fy(:,:,:,5) = TM*vM - KM*B2M
            Fy(:,:,:,6) = vM*rB1M - uM*rB2M
            Fy(:,:,:,7) = 0
            Fy(:,:,:,8) = vM*B3M - wM*B2M
                
            Fz(:,:,:,1) = uGint3D(:,:,:,4)
            Fz(:,:,:,2) = rhoM*uM*wM - B1M*B3M
            Fz(:,:,:,3) = rhoM*vM*wM - B2M*B3M
            Fz(:,:,:,4) = rhoM*wM**2 + SM - B3M**2
            Fz(:,:,:,5) = TM*wM - KM*B3M
            Fz(:,:,:,6) = wM*rB1M - uM*rB3M
            Fz(:,:,:,7) = wM*rB2M - vM*rB3M
            Fz(:,:,:,8) = 0
            
            
            RHS = 0
            RHS(:,:,:,2) = Fz(:,:,:,4)
            RHS(:,:,:,4) = -Fx(:,:,:,4)
            RHS(:,:,:,1:5) = RHS(:,:,:,1:5) - Fx(:,:,:,1:5)
            RHS(:,:,:,8) = 0
            
            RHSC = 0
            if (RHSCopen == 1) then
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        do kk = 0,Nphi
                            RHSC(i1,j1,kk,1) = S1(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                            RHSC(i1,j1,kk,2) = S2(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                            RHSC(i1,j1,kk,3) = S3(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                            RHSC(i1,j1,kk,4) = S4(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                            RHSC(i1,j1,kk,5) = S5(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK,gamma)
                            RHSC(i1,j1,kk,6) = (Rc(i) + hr1*lambda(i1))*S6(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                            RHSC(i1,j1,kk,7) = (Rc(i) + hr1*lambda(i1))*S7(Rc(i) + hr1*lambda(i1),Zc(j) + hz1*lambda(j1),Phi(kk),tRK)
                        end do
                    end do
                end do
            end if
                
            ! Fourier spectral method to calculate dH/dphi
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do n = 1,NumEq
                        Fzsin = 0
                        Fzcos = 0
                        Fzzsin = 0
                        Fzzcos = 0
                        do d = 1,Lphi
                            do kk = 0,Nphi
                                Fzsin(d) = Fzsin(d) + sink(kk,d)*Fz(i1,j1,kk,n)
                                Fzcos(d) = Fzcos(d) + cosk(kk,d)*Fz(i1,j1,kk,n)
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
                                RHS(i1,j1,kk,n) = RHS(i1,j1,kk,n) - Fzzsin(d)*sink(kk,d) - Fzzcos(d)*cosk(kk,d)
                            end do
                        end do
                        
                    end do
                end do
            end do
            
            RHS = RHS/(Rc(i) + hr1*RG) + RHSC
            
            do d = 1,dimPk1
                do n = 1,NumEq
                    do kk = 0,Nphi
                        do i1 = 1,NumGLP
                            do j1 = 1,NumGLP
                                if (d > 1) then
                                    du(i,j,kk,d,n) = du(i,j,kk,d,n) + 0.25d0*weight(i1)*weight(j1)*(Fx(i1,j1,kk,n)*phixG(i1,j1,d) + Fy(i1,j1,kk,n)*phiyG(i1,j1,d))
                                end if
                                du(i,j,kk,d,n) = du(i,j,kk,d,n) + 0.25d0*weight(i1)*weight(j1)*RHS(i1,j1,kk,n)*phiG(i1,j1,d)
                            end do
                        end do
                    end do
                end do
            end do
            
        end do
    end do
    !$omp end do
    
    ! The x-Flux
    !$omp single
    du8 = du(:,:,:,:,8)
    UR = 0
    UL = 0
    
    do i = 0,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        UR(i,j,kk,:,n) = UR(i,j,kk,:,n) + uh(i,j,kk,d,n)*phiGR(:,d)
                        UL(i + 1,j,kk,:,n) = UL(i + 1,j,kk,:,n) + uh(i + 1,j,kk,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do
    !$omp end single
    
    !$omp workshare
    rhoMR = uR(:,:,:,:,1)
    uMR = uR(:,:,:,:,2)/rhoMR
    vMR = uR(:,:,:,:,3)/rhoMR
    wMR = uR(:,:,:,:,4)/rhoMR
    EMR = uR(:,:,:,:,5)
    rB1MR = uR(:,:,:,:,6)
    rB2MR = uR(:,:,:,:,7)
    B1MR = rB1MR/Redge1(0:Nx,:,:,:)
    B2MR = rB2MR/Redge1(0:Nx,:,:,:)
    B3MR = uR(:,:,:,:,8)
    rB3MR = B3MR*Redge1(0:Nx,:,:,:)
    
    pMR = gamma1*(EMR - 0.5d0*rhoMR*(uMR**2 + vMR**2 + wMR**2) - 0.5d0*(B1MR**2 + B2MR**2 + B3MR**2))
    
    SMR = pMR + 0.5d0*(B1MR**2 + B2MR**2 + B3MR**2)
    TMR = EMR + SMR
    KMR = uMR*B1MR + vMR*B2MR + wMR*B3MR
    
    FR(:,:,:,:,1) = uR(:,:,:,:,2)
    FR(:,:,:,:,2) = rhoMR*uMR**2 + SMR - B1MR**2
    FR(:,:,:,:,3) = rhoMR*uMR*vMR - B1MR*B2MR
    FR(:,:,:,:,4) = rhoMR*uMR*wMR - B1MR*B3MR
    FR(:,:,:,:,5) = TMR*uMR - KMR*B1MR
    FR(:,:,:,:,6) = 0
    FR(:,:,:,:,7) = uMR*rB2MR - vMR*rB1MR
    FR(:,:,:,:,8) = uMR*B3MR - wMR*B1MR
    
    rhoML = uL(:,:,:,:,1)
    uML = uL(:,:,:,:,2)/rhoML
    vML = uL(:,:,:,:,3)/rhoML
    wML = uL(:,:,:,:,4)/rhoML
    EML = uL(:,:,:,:,5)
    rB1ML = uL(:,:,:,:,6)
    rB2ML = uL(:,:,:,:,7)
    B1ML = rB1ML/Redge1(0:Nx,:,:,:)
    B2ML = rB2ML/Redge1(0:Nx,:,:,:)
    B3ML = uL(:,:,:,:,8)
    rB3ML = B3ML*Redge1(0:Nx,:,:,:)
    
    pML = gamma1*(EML - 0.5d0*rhoML*(uML**2 + vML**2 + wML**2) - 0.5d0*(B1ML**2 + B2ML**2 + B3ML**2))
    
    SML = pML + 0.5d0*(B1ML**2 + B2ML**2 + B3ML**2)
    TML = EML + SML
    KML = uML*B1ML + vML*B2ML + wML*B3ML
    
    FL(:,:,:,:,1) = uL(:,:,:,:,2)
    FL(:,:,:,:,2) = rhoML*uML**2 + SML - B1ML**2
    FL(:,:,:,:,3) = rhoML*uML*vML - B1ML*B2ML
    FL(:,:,:,:,4) = rhoML*uML*wML - B1ML*B3ML
    FL(:,:,:,:,5) = TML*uML - KML*B1ML
    FL(:,:,:,:,6) = 0
    FL(:,:,:,:,7) = uML*rB2ML - vML*rB1ML
    FL(:,:,:,:,8) = uML*B3ML - wML*B1ML
    !$omp end workshare
    
    ! The y-Flux
    !$omp single
    UU = 0
    UD = 0
    
    do i = 1,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        UU(i,j,kk,:,n) = UU(i,j,kk,:,n) + uh(i,j,kk,d,n)*phiGU(:,d)
                        UD(i,j + 1,kk,:,n) = UD(i,j + 1,kk,:,n) + uh(i,j + 1,kk,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do
    !$omp end single
    
    !$omp workshare
    rhoMU = UU(:,:,:,:,1)
    uMU = UU(:,:,:,:,2)/rhoMU
    vMU = UU(:,:,:,:,3)/rhoMU
    wMU = UU(:,:,:,:,4)/rhoMU
    EMU = UU(:,:,:,:,5)
    rB1MU = UU(:,:,:,:,6)
    rB2MU = UU(:,:,:,:,7)
    B1MU = rB1MU/Redge2(:,0:Ny,:,:)
    B2MU = rB2MU/Redge2(:,0:Ny,:,:)
    B3MU = UU(:,:,:,:,8)
    rB3MU = B3MU*Redge2(:,0:Ny,:,:)
    
    pMU = gamma1*(EMU - 0.5d0*rhoMU*(uMU**2 + vMU**2 + wMU**2) - 0.5d0*(B1MU**2 + B2MU**2 + B3MU**2))
    
    SMU = pMU + 0.5d0*(B1MU**2 + B2MU**2 + B3MU**2)
    TMU = EMU + SMU
    KMU = uMU*B1MU + vMU*B2MU + wMU*B3MU
    
    FU(:,:,:,:,1) = uU(:,:,:,:,3)
    FU(:,:,:,:,2) = rhoMU*uMU*vMU - B1MU*B2MU
    FU(:,:,:,:,3) = rhoMU*vMU**2 + SMU - B2MU**2
    FU(:,:,:,:,4) = rhoMU*vMU*wMU - B2MU*B3MU
    FU(:,:,:,:,5) = TMU*vMU - KMU*B2MU
    FU(:,:,:,:,6) = vMU*rB1MU - uMU*rB2MU
    FU(:,:,:,:,7) = 0
    FU(:,:,:,:,8) = vMU*B3MU - wMU*B2MU
    
    rhoMD = UD(:,:,:,:,1)
    uMD = UD(:,:,:,:,2)/rhoMD
    vMD = UD(:,:,:,:,3)/rhoMD
    wMD = UD(:,:,:,:,4)/rhoMD
    EMD = UD(:,:,:,:,5)
    rB1MD = UD(:,:,:,:,6)
    rB2MD = UD(:,:,:,:,7)
    B1MD = rB1MD/Redge2(:,0:Ny,:,:)
    B2MD = rB2MD/Redge2(:,0:Ny,:,:)
    B3MD = UD(:,:,:,:,8)
    rB3MD = B3MD*Redge2(:,0:Ny,:,:)
    
    pMD = gamma1*(EMD - 0.5d0*rhoMD*(uMD**2 + vMD**2 - wMD**2) - 0.5d0*(B1MD**2 + B2MD**2 + B3MD**2))
    
    SMD = pMD + 0.5d0*(B1MD**2 + B2MD**2 + B3MD**2)
    TMD = EMD + SMD
    KMD = uMD*B1MD + vMD*B2MD + wMD*B3MD
    
    FD(:,:,:,:,1) = uD(:,:,:,:,3)
    FD(:,:,:,:,2) = rhoMD*uMD*vMD - B1MD*B2MD
    FD(:,:,:,:,3) = rhoMD*vMD**2 + SMD - B2MD**2
    FD(:,:,:,:,4) = rhoMD*vMD*wMD - B2MD*B3MD
    FD(:,:,:,:,5) = TMD*vMD - KMD*B2MD
    FD(:,:,:,:,6) = vMD*rB1MD - uMD*rB2MD
    FD(:,:,:,:,7) = 0
    FD(:,:,:,:,8) = vMD*B3MD - wMD*B2MD
    !$omp end workshare
    
    ! calculate Fx hat
    !$omp single
    do i = 0,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do j1 = 1,NumGLP
                    call eigenvalueMm(SRmax,SRmin,UR(i,j,kk,j1,1),UR(i,j,kk,j1,2),UR(i,j,kk,j1,3),UR(i,j,kk,j1,4),UR(i,j,kk,j1,5),B1MR(i,j,kk,j1),B2MR(i,j,kk,j1),UR(i,j,kk,j1,8),1,0)
                    call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,kk,j1,1),UL(i + 1,j,kk,j1,2),UL(i + 1,j,kk,j1,3),UL(i + 1,j,kk,j1,4),UL(i + 1,j,kk,j1,5),B1ML(i + 1,j,kk,j1),B2ML(i + 1,j,kk,j1),UL(i + 1,j,kk,j1,8),1,0)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FL(i + 1,j,kk,j1,:)
                    FL1 = FR(i,j,kk,j1,:)
                    UR1 = UL(i + 1,j,kk,j1,:)
                    UL1 = UR(i,j,kk,j1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 1
                        !call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 1
                        !call HLLD_Flux
                    else if (flux_type == 0) then
                        Fhat1 = 0.5d0*(FR1 + FL1)
                    end if
                    Fxhat(i,j,kk,j1,:) = Fhat1
                end do
            end do
        end do
    end do
    
    ! calculate Fy hat
    do i = 1,Nx
        do j = 0,Ny
            do kk = 0,Nphi
                do i1 = 1,NumGLP
                    call eigenvalueMm(SRmax,SRmin,UU(i,j,kk,i1,1),UU(i,j,kk,i1,2),UU(i,j,kk,i1,3),UU(i,j,kk,i1,4),UU(i,j,kk,i1,5),B1MU(i,j,kk,i1),B2MU(i,j,kk,i1),UU(i,j,kk,i1,8),0,1)
                    call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,kk,i1,1),UD(i,j + 1,kk,i1,2),UD(i,j + 1,kk,i1,3),UD(i,j + 1,kk,i1,4),UD(i,j + 1,kk,i1,5),B1MD(i,j + 1,kk,i1),B2MD(i,j + 1,kk,i1),UD(i,j + 1,kk,i1,8),0,1)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FD(i,j + 1,kk,i1,:)
                    FL1 = FU(i,j,kk,i1,:)
                    UR1 = UD(i,j + 1,kk,i1,:)
                    UL1 = UU(i,j,kk,i1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    else if (flux_type == 3) then
                        direction = 2
                        !call HLLC_Flux
                    else if (flux_type == 4) then
                        direction = 2
                        !call HLLD_Flux
                    else if (flux_type == 0) then
                        Fhat1 = 0.5d0*(FR1 + FL1)
                    end if
                    Fyhat(i,j,kk,i1,:) = Fhat1
                end do
            end do
        end do
    end do
    !$omp end single
    
    !!$omp single
    !$omp do collapse(6)
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        do j1 = 1,NumGLP
                            du(i,j,kk,d,n) = du(i,j,kk,d,n) - (0.5d0/hr)*weight(j1)*(Fxhat(i,j,kk,j1,n)*phiGR(j1,d) - Fxhat(i - 1,j,kk,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    
    !$omp do collapse(6)
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                do d = 1,dimPk1
                    do n = 1,NumEq
                        do i1 = 1,NumGLP
                            du(i,j,kk,d,n) = du(i,j,kk,d,n) - (0.5d0/hz)*weight(i1)*(Fyhat(i,j,kk,i1,n)*phiGU(i1,d) - Fyhat(i,j - 1,kk,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    
    !$omp end parallel
    
    do d = 1,dimPk1
        du(:,:,:,d,:) = du(:,:,:,d,:)/mm(d)
    end do
    
    end subroutine Lh