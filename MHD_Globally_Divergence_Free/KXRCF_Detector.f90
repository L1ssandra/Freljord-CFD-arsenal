    subroutine KXRCF_Detector
    
    include 'com.txt'
    
    real Iij(NumEq),Lij,uinfty(NumEq),Ck
    
    Ck = 0
    
    !$omp parallel num_threads(12) default(shared)
    
    uG = 0
    !$omp do collapse(4)
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,NumEq
                do d = 1,dimPk
                    uG(i,j,:,:,n) = uG(i,j,:,:,n) + uh(i,j,d,n)*phiG(:,:,d)
                end do
            end do
        end do
    end do
    !$omp end do
    
    UR = 0
    UL = 0
    !$omp do collapse(4)
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UR(i,j,:,n) = UR(i,j,:,n) + uh(i,j,d,n)*phiGR(:,d)
                    UL(i + 1,j,:,n) = UL(i + 1,j,:,n) + uh(i + 1,j,d,n)*phiGL(:,d)
                end do
            end do
        end do
    end do
    !$omp end do
    
    UU = 0
    UD = 0
    !$omp do collapse(4)
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,dimPk
                do n = 1,NumEq
                    UU(i,j,:,n) = UU(i,j,:,n) + uh(i,j,d,n)*phiGU(:,d)
                    UD(i,j + 1,:,n) = UD(i,j + 1,:,n) + uh(i,j + 1,d,n)*phiGD(:,d)
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
    
    do i = 1,Nx
        do j = 1,Ny
            Iij = 0
            Lij = 1e-12
            if (uh(i + 1,j,1,2) < 0) then ! Rightside inflow
                do j1 = 1,NumGLP
                    Iij = Iij + hy1*weight(j1)*(UR(i,j,j1,:) - UL(i + 1,j,j1,:))
                    Lij = Lij + hy
                end do
            end if
            
            if (uh(i - 1,j,1,2) > 0) then ! Leftside inflow
                do j1 = 1,NumGLP
                    Iij = Iij + hy1*weight(j1)*(UL(i,j,j1,:) - UR(i - 1,j,j1,:))
                    Lij = Lij + hy
                end do
            end if
            
            if (uh(i,j + 1,1,3) < 0) then ! Upside inflow
                do i1 = 1,NumGLP
                    Iij = Iij + hx1*weight(i1)*(UU(i,j,i1,:) - UD(i,j + 1,i1,:))
                    Lij = Lij + hx
                end do
            end if
            
            if (uh(i,j - 1,1,3) > 0) then ! Downside inflow
                do i1 = 1,NumGLP
                    Iij = Iij + hx1*weight(i1)*(UD(i,j,i1,:) - UU(i,j - 1,i1,:))
                    Lij = Lij + hx
                end do
            end if
            
            uinfty = 0
            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    do n = 1,NumEq
                        if (abs(uG(i,j,i1,j1,n)) + 1e-12 > uinfty(n)) then
                            uinfty(n) = abs(uG(i,j,i1,j1,n)) + 1e-12
                        end if
                    end do
                end do
            end do
            
            Iij = abs(Iij)/((hx*hy)**(0.25*(k + 1))*Lij*uinfty)
            
            do n = 1,NumEq
                if (Iij(n) >= Ck) then
                    Is_Trouble_Cell(i,j) = 1
                end if
            end do
            
        end do
    end do

    end subroutine KXRCF_Detector
            
                