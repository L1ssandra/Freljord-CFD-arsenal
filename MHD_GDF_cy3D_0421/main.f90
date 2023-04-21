    program main

    include 'com.txt'
    
    call get_GLP
    
    call init_data
    
    call RK3
    
    call calculate_L2_error
    
    call save_RZ
    
    call save_solution

    end program main

