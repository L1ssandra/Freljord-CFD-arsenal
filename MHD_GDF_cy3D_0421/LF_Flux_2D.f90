    subroutine LF_Flux_2D
    
    include 'com.txt'
    
    real alphaxRU,alphaxLU,alphaxRD,alphaxLD
    real alphayRU,alphayLU,alphayRD,alphayLD
    real alphax2D,alphay2D
    real EzRU,EzLU,EzRD,EzLD
    real B1RU,B1LU,B1RD,B1LD
    real B2RU,B2LU,B2RD,B2LD
    real EzR1,EzL1,EzU1,EzD1
    
    alphax2D = max(abs(URU1(2)/URU1(1)),abs(ULU1(2)/ULU1(1)),abs(URD1(2)/URD1(1)),abs(ULD1(2)/ULD1(1)))
    alphay2D = max(abs(URU1(3)/URU1(1)),abs(ULU1(3)/ULU1(1)),abs(URD1(3)/URD1(1)),abs(ULD1(3)/ULD1(1)))
    
    call calculate_Ez(EzRU,URU1(2)/URU1(1),URU1(3)/URU1(1),URU1(6),URU1(7))
    call calculate_Ez(EzLU,ULU1(2)/ULU1(1),ULU1(3)/ULU1(1),ULU1(6),ULU1(7))
    call calculate_Ez(EzRD,URD1(2)/URD1(1),URD1(3)/URD1(1),URD1(6),URD1(7))
    call calculate_Ez(EzLD,ULD1(2)/ULD1(1),ULD1(3)/ULD1(1),ULD1(6),ULD1(7))
    
    B1RU = URU1(6)
    B1LU = ULU1(6)
    B1RD = URD1(6)
    B1LD = ULD1(6)
    
    B2RU = URU1(7)
    B2LU = ULU1(7)
    B2RD = URD1(7)
    B2LD = ULD1(7)
    
    Ezhat = 0.25*(EzRU + EzLU + EzRD + EzLD) - 0.25*alphay2D*(0.5*(B1RU + B1LU) - 0.5*(B1RD + B1LD)) + 0.25*alphax2D*(0.5*(B2RU + B2RD) - 0.5*(B2LU + B2LD))
    
    end subroutine LF_Flux_2D