    function pS(r,z,phi)
    
    real r,z,phi,rr
    real pS
	
    rr = ((r - 0.5)**2 + z**2)**0.5
    
	if (rr > 0.1) then
        pS = 0.1
    else
        pS = 10
    end if
    
	end