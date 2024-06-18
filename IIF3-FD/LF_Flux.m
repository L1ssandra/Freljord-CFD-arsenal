function fhat = LF_Flux(uR,uL,fR,fL,alpha)

fhat = 0.5*(fR + fL - alpha*(uR - uL));

end