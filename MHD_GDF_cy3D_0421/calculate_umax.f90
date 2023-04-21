    subroutine calculate_umax
    
    include 'com.txt'
    
    umax = 0
    
    do i = 1,Nx
        do j = 1,Ny
            do kk = 0,Nphi
                if (abs(uh(i,j,kk,1,1)) > umax) then
                    umax = abs(uh(i,j,kk,1,1))
                end if
            end do
        end do
    end do
    
    end subroutine calculate_umax