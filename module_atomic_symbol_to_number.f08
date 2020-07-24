 MODULE module_atomic_symbol_to_number
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 IMPLICIT NONE
 
  CONTAINS
 
 SUBROUTINE atomic_symbol_to_number(atomic_symbol,output_FID,z)
 !-----------------------------------------------------------------------------------
 
 INTEGER :: z,output_FID
 CHARACTER(2) :: atomic_symbol
 
 IF (trim(atomic_symbol) == 'H') THEN
    z = 1
 ELSE IF (trim(atomic_symbol) == 'He') THEN
    z = 2
 ELSE IF (trim(atomic_symbol) == 'Li') THEN
    z = 3
 ELSE IF (trim(atomic_symbol) == 'Be') THEN
    z = 4
 ELSE IF (trim(atomic_symbol) == 'B') THEN
    z = 5
 ELSE IF (trim(atomic_symbol) == 'C') THEN
    z = 6
 ELSE IF (trim(atomic_symbol) == 'N') THEN
    z = 7
 ELSE IF (trim(atomic_symbol) == 'O') THEN
    z = 8
 ELSE IF (trim(atomic_symbol) == 'F') THEN
    z = 9
 ELSE IF (trim(atomic_symbol) == 'Ne') THEN
    z = 10
 ELSE IF (trim(atomic_symbol) == 'Na') THEN
    z = 11 
 ELSE IF (trim(atomic_symbol) == 'Mg') THEN
    z = 12 
 ELSE IF (trim(atomic_symbol) == 'Al') THEN
    z = 13 
 ELSE IF (trim(atomic_symbol) == 'Si') THEN
    z = 14 
 ELSE IF (trim(atomic_symbol) == 'P') THEN
    z = 15 
 ELSE IF (trim(atomic_symbol) == 'S') THEN
    z = 16 
 ELSE IF (trim(atomic_symbol) == 'Cl') THEN
    z = 17 
 ELSE IF (trim(atomic_symbol) == 'Ar') THEN
    z = 18 
 ELSE IF (trim(atomic_symbol) == 'K') THEN
    z = 19 
 ELSE IF (trim(atomic_symbol) == 'Ca') THEN
    z = 20 
 ELSE IF (trim(atomic_symbol) == 'Sc') THEN
    z = 21 
 ELSE IF (trim(atomic_symbol) == 'Ti') THEN
    z = 22 
 ELSE IF (trim(atomic_symbol) == 'V') THEN
    z = 23 
 ELSE IF (trim(atomic_symbol) == 'Cr') THEN
    z = 24 
 ELSE IF (trim(atomic_symbol) == 'Mn') THEN
    z = 25 
 ELSE IF (trim(atomic_symbol) == 'Fe') THEN
    z = 26 
 ELSE IF (trim(atomic_symbol) == 'Co') THEN
    z = 27 
 ELSE IF (trim(atomic_symbol) == 'Ni') THEN
    z = 28 
 ELSE IF (trim(atomic_symbol) == 'Cu') THEN
    z = 29 
 ELSE IF (trim(atomic_symbol) == 'Zn') THEN
    z = 30 
 ELSE IF (trim(atomic_symbol) == 'Ga') THEN
    z = 31 
 ELSE IF (trim(atomic_symbol) == 'Ge') THEN
    z = 32 
 ELSE IF (trim(atomic_symbol) == 'As') THEN
    z = 33 
 ELSE IF (trim(atomic_symbol) == 'Se') THEN
    z = 34 
 ELSE IF (trim(atomic_symbol) == 'Br') THEN
    z = 35 
 ELSE IF (trim(atomic_symbol) == 'Kr') THEN
    z = 36 
 ELSE IF (trim(atomic_symbol) == 'Rb') THEN
    z = 37 
 ELSE IF (trim(atomic_symbol) == 'Sr') THEN
    z = 38 
 ELSE IF (trim(atomic_symbol) == 'Y') THEN
    z = 39 
 ELSE IF (trim(atomic_symbol) == 'Zr') THEN
    z = 40 
 ELSE IF (trim(atomic_symbol) == 'Nb') THEN
    z = 41 
 ELSE IF (trim(atomic_symbol) == 'Mo') THEN
    z = 42 
 ELSE IF (trim(atomic_symbol) == 'Tc') THEN
    z = 43 
 ELSE IF (trim(atomic_symbol) == 'Ru') THEN
    z = 44 
 ELSE IF (trim(atomic_symbol) == 'Rh') THEN
    z = 45 
 ELSE IF (trim(atomic_symbol) == 'Pd') THEN
    z = 46 
 ELSE IF (trim(atomic_symbol) == 'Ag') THEN
    z = 47 
 ELSE IF (trim(atomic_symbol) == 'Cd') THEN
    z = 48 
 ELSE IF (trim(atomic_symbol) == 'In') THEN
    z = 49 
 ELSE IF (trim(atomic_symbol) == 'Sn') THEN
    z = 50 
 ELSE IF (trim(atomic_symbol) == 'Sb') THEN
    z = 51 
 ELSE IF (trim(atomic_symbol) == 'Te') THEN
    z = 52 
 ELSE IF (trim(atomic_symbol) == 'I') THEN
    z = 53 
 ELSE IF (trim(atomic_symbol) == 'Xe') THEN
    z = 54 
 ELSE IF (trim(atomic_symbol) == 'Cs') THEN
    z = 55 
 ELSE IF (trim(atomic_symbol) == 'Ba') THEN
    z = 56 
 ELSE IF (trim(atomic_symbol) == 'La') THEN
    z = 57 
 ELSE IF (trim(atomic_symbol) == 'Ce') THEN
    z = 58 
 ELSE IF (trim(atomic_symbol) == 'Pr') THEN
    z = 59 
 ELSE IF (trim(atomic_symbol) == 'Nd') THEN
    z = 60 
 ELSE IF (trim(atomic_symbol) == 'Pm') THEN
    z = 61 
 ELSE IF (trim(atomic_symbol) == 'Sm') THEN
    z = 62 
 ELSE IF (trim(atomic_symbol) == 'Eu') THEN
    z = 63 
 ELSE IF (trim(atomic_symbol) == 'Gd') THEN
    z = 64 
 ELSE IF (trim(atomic_symbol) == 'Tb') THEN
    z = 65 
 ELSE IF (trim(atomic_symbol) == 'Dy') THEN
    z = 66 
 ELSE IF (trim(atomic_symbol) == 'Ho') THEN
    z = 67 
 ELSE IF (trim(atomic_symbol) == 'Er') THEN
    z = 68 
 ELSE IF (trim(atomic_symbol) == 'Tm') THEN
    z = 69 
 ELSE IF (trim(atomic_symbol) == 'Yb') THEN
    z = 70 
 ELSE IF (trim(atomic_symbol) == 'Lu') THEN
    z = 71 
 ELSE IF (trim(atomic_symbol) == 'Hf') THEN
    z = 72 
 ELSE IF (trim(atomic_symbol) == 'Ta') THEN
    z = 73 
 ELSE IF (trim(atomic_symbol) == 'W') THEN
    z = 74 
 ELSE IF (trim(atomic_symbol) == 'Re') THEN
    z = 75 
 ELSE IF (trim(atomic_symbol) == 'Os') THEN
    z = 76 
 ELSE IF (trim(atomic_symbol) == 'Ir') THEN
    z = 77 
 ELSE IF (trim(atomic_symbol) == 'Pt') THEN
    z = 78 
 ELSE IF (trim(atomic_symbol) == 'Au') THEN
    z = 79 
 ELSE IF (trim(atomic_symbol) == 'Hg') THEN
    z = 80 
 ELSE IF (trim(atomic_symbol) == 'Tl') THEN
    z = 81 
 ELSE IF (trim(atomic_symbol) == 'Pb') THEN
    z = 82 
 ELSE IF (trim(atomic_symbol) == 'Bi') THEN
    z = 83 
 ELSE IF (trim(atomic_symbol) == 'Po') THEN
    z = 84 
 ELSE IF (trim(atomic_symbol) == 'At') THEN
    z = 85 
 ELSE IF (trim(atomic_symbol) == 'Rn') THEN
    z = 86 
 ELSE IF (trim(atomic_symbol) == 'Fr') THEN
    z = 87 
 ELSE IF (trim(atomic_symbol) == 'Ra') THEN
    z = 88 
 ELSE IF (trim(atomic_symbol) == 'Ac') THEN
    z = 89 
 ELSE IF (trim(atomic_symbol) == 'Th') THEN
    z = 90 
 ELSE IF (trim(atomic_symbol) == 'Pa') THEN
    z = 91 
 ELSE IF (trim(atomic_symbol) == 'U') THEN
    z = 92 
 ELSE IF (trim(atomic_symbol) == 'Np') THEN
    z = 93 
 ELSE IF (trim(atomic_symbol) == 'Pu') THEN
    z = 94 
 ELSE IF (trim(atomic_symbol) == 'Am') THEN
    z = 95 
 ELSE IF (trim(atomic_symbol) == 'Cm') THEN
    z = 96 
 ELSE IF (trim(atomic_symbol) == 'Bk') THEN
    z = 97 
 ELSE IF (trim(atomic_symbol) == 'Cf') THEN
    z = 98 
 ELSE IF (trim(atomic_symbol) == 'Es') THEN
    z = 99 
 ELSE IF (trim(atomic_symbol) == 'Fm') THEN
    z = 100 
 ELSE IF (trim(atomic_symbol) == 'Md') THEN
    z = 101 
 ELSE IF (trim(atomic_symbol) == 'No') THEN
    z = 102 
 ELSE IF (trim(atomic_symbol) == 'Lr') THEN
    z = 103 
 ELSE IF (trim(atomic_symbol) == 'Rf') THEN
    z = 104 
 ELSE IF (trim(atomic_symbol) == 'Db') THEN
    z = 105 
 ELSE IF (trim(atomic_symbol) == 'Sg') THEN
    z = 106 
 ELSE IF (trim(trim(atomic_symbol)) == 'Bh') THEN
    z = 107 
 ELSE IF (trim(trim(atomic_symbol)) == 'Hs') THEN
    z = 108 
 ELSE IF (trim(trim(atomic_symbol)) == 'Mt') THEN
    z = 109 
 ! Alternate spellings
 ELSE IF (trim(trim(atomic_symbol)) == 'X') THEN
    z = 54 
 ELSE
    WRITE(output_FID,*)'Atomic symbol not properly identified. Program will terminate.'
    STOP
 END IF

 END SUBROUTINE atomic_symbol_to_number

 END MODULE module_atomic_symbol_to_number