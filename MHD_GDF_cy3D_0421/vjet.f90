    function vjet(r,z,phi)
    
    real r,z,phi
    real vjet
    
	if (r <= 1.5 .and. z <= 3) then
        vjet = 1
    else
        vjet = 0
    end if
    
    end
    
    
    function rhojet(r,z,phi)
    
    real r,z,phi
    real rhojet
    
	if (r <= 1.5 .and. z <= 3) then
        rhojet = 100
    else
        rhojet = 10
    end if
    
    end