    subroutine calculate_dt
    
    include 'com.txt'
    
    if (RKorder == 1) then
        CFL = 0.05
    else if (RKorder == 3) then
        CFL = 0.2
    end if
    
    alphax = 0
    alphay = 0
    
    do i = 1,Nr
        do j = 1,Nz
            do kk = 0,Nphi
                call eigenvalueMm(alpha1,alpha2,uh(i,j,kk,1,1),uh(i,j,kk,1,2),uh(i,j,kk,1,3),uh(i,j,kk,1,4),uh(i,j,kk,1,5),uh(i,j,kk,1,6)/Rc(i),uh(i,j,kk,1,7)/Rc(i),uh(i,j,kk,1,8),1,0)
                if (abs(alpha1) > alphax .or. abs(alpha2) > alphax) then
                    alphax = max(abs(alpha1),abs(alpha2))
                end if
                call eigenvalueMm(alpha1,alpha2,uh(i,j,kk,1,1),uh(i,j,kk,1,2),uh(i,j,kk,1,3),uh(i,j,kk,1,4),uh(i,j,kk,1,5),uh(i,j,kk,1,6)/Rc(i),uh(i,j,kk,1,7)/Rc(i),uh(i,j,kk,1,8),0,1)
                if (abs(alpha1) > alphay .or. abs(alpha2) > alphay) then
                    alphay = max(abs(alpha1),abs(alpha2))
                end if
            end do
        end do
    end do
    
    dt = CFL/(alphax/hr + alphay/hz)
    
    end subroutine calculate_dt