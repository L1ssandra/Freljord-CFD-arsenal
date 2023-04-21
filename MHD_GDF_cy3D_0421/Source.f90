    function S1(r,z,phi,t)
    
    real r,z,phi,t
    real S1
	
    S1 = (1 + 0.5*sin(r + z - 2*t))/r
    
	end
	!*************************************************
    function S2(r,z,phi,t)
    
    real r,z,phi,t
    real S2
	
    S2 = 0.25*sin(r - t)*cos(r - t)*cos(phi)**2 + 0.25*sin(r - t)*cos(z - t)*cos(phi)**2 + (1 + 0.5*sin(r + z - 2*t))/r
    
	end
	!*************************************************
	function S3(r,z,phi,t)
    
    real r,z,phi,t
    real S3
	
    S3 = 0.25*cos(r - t)*sin(z - t)*cos(phi)**2 + 0.25*sin(z - t)*cos(z - t) + (1 + 0.5*sin(r + z - 2*t))/r + 0.25*sin(r - t)*sin(z - t)*sin(phi)**2/r
    
	end
	!*************************************************
	function S4(r,z,phi,t)
    
    real r,z,phi,t
    real S4
	
    S4 = -0.25*sin(r - t)*cos(z - t)*sin(phi)*cos(phi) - 0.25*sin(r - t)**2*sin(phi)*cos(phi)/r
    
	end
	!*************************************************
	function S5(r,z,phi,t,gamma)
    
    real r,z,phi,t
    real S5
	
    S5 = 0.25*sin(r - t)*cos(r - t)*cos(phi)**2 + 0.25*cos(r - t)*sin(z - t)*cos(phi)**2 + 0.25*sin(z - t)*cos(z - t) + 0.25*sin(r - t)*cos(z - t)*cos(phi)**2 + 0.25*(sin(r - t) - sin(z - t))*sin(z - t)*sin(phi)**2/r + (gamma/(gamma - 1) + 1 + 0.5*sin(r + z - 2*t) + 0.25*sin(z - t)**2 + 0.25*sin(r - t)**2*cos(phi)**2)/r
    
    end
    !*************************************************
	function S6(r,z,phi,t)
    
    real r,z,phi,t
    real S6
	
    S6 = -0.5*sin(z - t)*cos(phi)/r
    
    end
    !*************************************************
	function S7(r,z,phi,t)
    
    real r,z,phi,t
    real S7
	
    S7 = 0.5*sin(r - t)*cos(phi)/r
    
	end