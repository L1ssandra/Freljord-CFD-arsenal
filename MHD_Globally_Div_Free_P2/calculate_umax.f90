    subroutine calculate_umax
    
    include 'com.txt'
    
    umax = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do n = 1,numEq
                if (abs(uh(i,j,1,n)) > umax) then
                    umax = abs(uh(i,j,1,n))
                end if
            end do
        end do
    end do
    
    end subroutine calculate_umax