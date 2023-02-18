    subroutine L2_Pro
    
    include 'com.txt'
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,numEq
                uh(i,j,:,n) = 0
                do d = 1,dimPk
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            if (d <= 6 .or. (n == 6 .and. (d == 7 .or. d == 9)) .or. (n == 7 .and. (d == 8 .or. d == 10))) then
                                uh(i,j,d,n) = uh(i,j,d,n) + weight(i1)*weight(j1)*phiG(i1,j1,d)*ureal(i,j,i1,j1,n)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    
    uh = 0.25*uh
    
    do d = 1,dimPk
        uh(:,:,d,:) = uh(:,:,d,:)/mm(d)
    end do
    
    if (shock == 1) then
        do i = 1,Nx
            do j = 1,Ny
                uh(i,j,2:dimPk,:) = 0
                uh(i,j,1,:) = ureal(i,j,(NumGLP + 1)/2,(NumGLP + 1)/2,:)
            end do
        end do
    end if
    
    Bx = 0
    By = 0
    do i = 0,Nx
        do j = 1,Ny
            do d = 1,k + 1
                do j1 = 1,NumGLP
                    Bx(i,j,d) = Bx(i,j,d) + 0.5*weight(j1)*EzG(j1,d)*Ez1real(i,j,j1)
                end do
            end do
        end do
    end do
    
    do i = 1,Nx
        do j = 0,Ny
            do d = 1,k + 1
                do i1 = 1,NumGLP
                    By(i,j,d) = By(i,j,d) + 0.5*weight(i1)*EzG(i1,d)*Ez2real(i,j,i1)
                end do
            end do
        end do
    end do
    
    do d = 1,k + 1
        Bx(:,:,d) = Bx(:,:,d)/mmE(d)
        By(:,:,d) = By(:,:,d)/mmE(d)
    end do
    
    end subroutine L2_Pro