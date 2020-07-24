 MODULE module_print_overlap_populations
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations
 USE module_string_utilities
 USE module_atomic_number_to_symbol

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE print_overlap_populations ()
 !===================================================================================
 ! Print the overlap populations of each atom pair
 !=================================================================================== 
 
 CHARACTER(29) :: filename 
 INTEGER :: zero_FID,j
 REAL(kind=dp) :: V1(3),V2(3),V3(3)
 
 !Generate the output file
 filename = 'overlap_populations.xyz'
 OPEN(NEWUNIT=zero_FID, FILE=filename,STATUS='REPLACE')
 WRITE(zero_FID,'(i5)') 2*num_sym_unique_bond_pairs
 FLUSH(zero_FID)
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 IF ((periodicA .or. periodicB) .or. periodicC) THEN
   WRITE(zero_FID,'(a,3f13.6,a,3f13.6,a,3f13.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell [{ ',V1(1),V1(2),&
   V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'    
   FLUSH(zero_FID)
 ELSE
   WRITE(zero_FID,'(a)') 'Nonperiodic system'
   FLUSH(zero_FID)
 END IF
 WRITE(zero_FID,'(a)') 'atom1, atom2, translation A, translation B, translation C, overlap population'
 DO j = 1,num_sym_unique_bond_pairs
   WRITE(zero_FID,'(i6,a,i6,a,i6,a,i6,a,i6,a,f17.10)') nint(bond_pair_matrix(1,j)),' ',nint(bond_pair_matrix(2,j)),' ',&
   nint(bond_pair_matrix(3,j)),' ',nint(bond_pair_matrix(4,j)),' ',nint(bond_pair_matrix(5,j)),' ',bond_pair_matrix(14,j)
   WRITE(zero_FID,'(i6,a,i6,a,i6,a,i6,a,i6,a,f17.10)') nint(bond_pair_matrix(2,j)),' ',nint(bond_pair_matrix(1,j)),' ',&
   -nint(bond_pair_matrix(3,j)),' ',-nint(bond_pair_matrix(4,j)),' ',-nint(bond_pair_matrix(5,j)),' ',bond_pair_matrix(14,j)
   FLUSH(zero_FID)
 END DO
 WRITE(zero_FID,*)' '
 CLOSE(zero_FID)
  
 END SUBROUTINE print_overlap_populations
 END MODULE module_print_overlap_populations