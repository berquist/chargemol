 MODULE module_atomic_number_to_symbol
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE atomic_number_to_symbol(z,output_fid,atomic_symbol)
 !=====================================================================================
 
 INTEGER :: z,output_FID
 CHARACTER(2) :: atomic_symbol
 atomic_symbol=''

 IF (z == 1) THEN
    atomic_symbol = 'H'
 ELSE IF (z == 2) THEN
    atomic_symbol = 'He'    
 ELSE IF (z == 3) THEN
    atomic_symbol = 'Li'
 ELSE IF (z == 4) THEN
    atomic_symbol = 'Be' 
 ELSE IF (z == 5) THEN
    atomic_symbol = 'B'
 ELSE IF (z == 6) THEN
    atomic_symbol = 'C'
 ELSE IF (z == 7) THEN
    atomic_symbol = 'N'
 ELSE IF (z == 8) THEN
    atomic_symbol = 'O'
 ELSE IF (z == 9) THEN
    atomic_symbol = 'F'
 ELSE IF (z == 10) THEN
    atomic_symbol = 'Ne'
 ELSE IF (z == 11) THEN
    atomic_symbol = 'Na'
 ELSE IF (z == 12) THEN
    atomic_symbol = 'Mg'
 ELSE IF (z == 13) THEN
    atomic_symbol = 'Al'
 ELSE IF (z == 14) THEN
    atomic_symbol = 'Si'
 ELSE IF (z == 15) THEN
    atomic_symbol = 'P'
 ELSE IF (z == 16) THEN
    atomic_symbol = 'S'
 ELSE IF (z == 17) THEN
    atomic_symbol = 'Cl'
 ELSE IF (z == 18) THEN
    atomic_symbol = 'Ar'
 ELSE IF (z == 19) THEN
    atomic_symbol = 'K'
 ELSE IF (z == 20) THEN
    atomic_symbol = 'Ca'
 ELSE IF (z == 21) THEN
    atomic_symbol = 'Sc'
 ELSE IF (z == 22) THEN
    atomic_symbol = 'Ti'
 ELSE IF (z == 23) THEN
    atomic_symbol = 'V'
 ELSE IF (z == 24) THEN
    atomic_symbol = 'Cr'
 ELSE IF (z == 25) THEN
    atomic_symbol = 'Mn'
 ELSE IF (z == 26) THEN
    atomic_symbol = 'Fe'
 ELSE IF (z == 27) THEN
    atomic_symbol = 'Co'
 ELSE IF (z == 28) THEN
    atomic_symbol = 'Ni'
 ELSE IF (z == 29) THEN
    atomic_symbol = 'Cu'
 ELSE IF (z == 30) THEN
    atomic_symbol = 'Zn'
 ELSE IF (z == 31) THEN
    atomic_symbol = 'Ga'
 ELSE IF (z == 32) THEN
    atomic_symbol = 'Ge'
 ELSE IF (z == 33) THEN
    atomic_symbol = 'As'
 ELSE IF (z == 34) THEN
    atomic_symbol = 'Se'
 ELSE IF (z == 35) THEN
    atomic_symbol = 'Br'
 ELSE IF (z == 36) THEN
    atomic_symbol = 'Kr'
 ELSE IF (z == 37) THEN
    atomic_symbol = 'Rb'
 ELSE IF (z == 38) THEN
    atomic_symbol = 'Sr'
 ELSE IF (z == 39) THEN
    atomic_symbol = 'Y'
 ELSE IF (z == 40) THEN
    atomic_symbol = 'Zr'
 ELSE IF (z == 41) THEN
    atomic_symbol = 'Nb'
 ELSE IF (z == 42) THEN
    atomic_symbol = 'Mo'
 ELSE IF (z == 43) THEN
    atomic_symbol = 'Tc'
 ELSE IF (z == 44) THEN
    atomic_symbol = 'Ru'
 ELSE IF (z == 45) THEN
    atomic_symbol = 'Rh'
 ELSE IF (z == 46) THEN
    atomic_symbol = 'Pd'
 ELSE IF (z == 47) THEN
    atomic_symbol = 'Ag'
 ELSE IF (z == 48) THEN
    atomic_symbol = 'Cd'
 ELSE IF (z == 49) THEN
    atomic_symbol = 'In'
 ELSE IF (z == 50) THEN
    atomic_symbol = 'Sn'
 ELSE IF (z == 51) THEN
    atomic_symbol = 'Sb'
 ELSE IF (z == 52) THEN
    atomic_symbol = 'Te'
 ELSE IF (z == 53) THEN 
    atomic_symbol = 'I'
 ELSE IF (z == 54) THEN
    atomic_symbol = 'Xe'
 ELSE IF (z == 55) THEN
    atomic_symbol = 'Cs'
 ELSE IF (z == 56) THEN
    atomic_symbol = 'Ba'
 ELSE IF (z == 57) THEN
    atomic_symbol = 'La'
 ELSE IF (z == 58) THEN
    atomic_symbol = 'Ce'
 ELSE IF (z == 59) THEN
    atomic_symbol = 'Pr'
 ELSE IF (z == 60) THEN
    atomic_symbol = 'Nd'
 ELSE IF (z == 61) THEN
    atomic_symbol = 'Pm'
 ELSE IF (z == 62) THEN
    atomic_symbol = 'Sm'
 ELSE IF (z == 63) THEN
    atomic_symbol = 'Eu'
 ELSE IF (z == 64) THEN
    atomic_symbol = 'Gd'
 ELSE IF (z == 65) THEN
    atomic_symbol = 'Tb'
 ELSE IF (z == 66) THEN
    atomic_symbol = 'Dy'
 ELSE IF (z == 67) THEN
    atomic_symbol = 'Ho'
 ELSE IF (z == 68) THEN
    atomic_symbol = 'Er'
 ELSE IF (z == 69) THEN
    atomic_symbol = 'Tm'
 ELSE IF (z == 70) THEN
    atomic_symbol = 'Yb'
 ELSE IF (z == 71) THEN
    atomic_symbol = 'Lu'
 ELSE IF (z == 72) THEN
    atomic_symbol = 'Hf'
 ELSE IF (z == 73) THEN
    atomic_symbol = 'Ta'
 ELSE IF (z == 74) THEN
    atomic_symbol = 'W'
 ELSE IF (z == 75) THEN
    atomic_symbol = 'Re'
 ELSE IF (z == 76) THEN
    atomic_symbol = 'Os'
 ELSE IF (z == 77) THEN
    atomic_symbol = 'Ir'
 ELSE IF (z == 78) THEN
    atomic_symbol = 'Pt'
 ELSE IF (z == 79) THEN
    atomic_symbol = 'Au'
 ELSE IF (z == 80) THEN
    atomic_symbol = 'Hg'
 ELSE IF (z == 81) THEN
    atomic_symbol = 'Tl'
 ELSE IF (z == 82) THEN
    atomic_symbol = 'Pb'
 ELSE IF (z == 83) THEN
    atomic_symbol = 'Bi'
 ELSE IF (z == 84) THEN
    atomic_symbol = 'Po'
 ELSE IF (z == 85) THEN
    atomic_symbol = 'At'
 ELSE IF (z == 86) THEN
    atomic_symbol = 'Rn'
 ELSE IF (z == 87) THEN
    atomic_symbol = 'Fr'
 ELSE IF (z == 88) THEN
    atomic_symbol = 'Ra'
 ELSE IF (z == 89) THEN
    atomic_symbol = 'Ac'
 ELSE IF (z == 90) THEN
    atomic_symbol = 'Th'
 ELSE IF (z == 91) THEN
    atomic_symbol = 'Pa'
 ELSE IF (z == 92) THEN
    atomic_symbol = 'U'
 ELSE IF (z == 93) THEN
    atomic_symbol = 'Np'
 ELSE IF (z == 94) THEN
    atomic_symbol = 'Pu'
 ELSE IF (z == 95) THEN
    atomic_symbol = 'Am'
 ELSE IF (z == 96) THEN
    atomic_symbol = 'Cm'
 ELSE IF (z == 97) THEN
    atomic_symbol = 'Bk'
 ELSE IF (z == 98) THEN
    atomic_symbol = 'Cf'
 ELSE IF (z == 99) THEN
    atomic_symbol = 'Es'
 ELSE IF (z == 100) THEN
    atomic_symbol = 'Fm'
 ELSE IF (z == 101) THEN
    atomic_symbol = 'Md'
 ELSE IF (z == 102) THEN
    atomic_symbol = 'No'
 ELSE IF (z == 103) THEN
    atomic_symbol = 'Lr'
 ELSE IF (z == 104) THEN
    atomic_symbol = 'Rf'
 ELSE IF (z == 105) THEN
    atomic_symbol = 'Db'
 ELSE IF (z == 106) THEN
    atomic_symbol = 'Sg'
 ELSE IF (z == 107) THEN
    atomic_symbol = 'Bh'
 ELSE IF (z == 108) THEN
    atomic_symbol = 'Hs'
 ELSE IF (z == 109) THEN
    atomic_symbol = 'Mt'
 ELSE
    WRITE(output_FID,*)'Atomic symbol not properly identified. Program will terminate.'
    FLUSH(output_FID)
    STOP
 END IF
 
 END SUBROUTINE atomic_number_to_symbol
 
 END MODULE module_atomic_number_to_symbol