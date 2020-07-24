 MODULE module_compute_center_of_mass
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_global_parameters
 USE module_common_variable_declarations

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE compute_center_of_mass()
 !=================================================================================== 
 
 REAL(kind=dp) :: sumvector(3),Mtot
 
 !standard_atomic_weights
 sumvector=0.0_dp
 Mtot=0.0_dp
 DO j=1,natoms
   sumvector(1) = sumvector(1) + atomic_weight(atomic_number(j))*coords(1,j)
   sumvector(2) = sumvector(2) + atomic_weight(atomic_number(j))*coords(2,j)
   sumvector(3) = sumvector(3) + atomic_weight(atomic_number(j))*coords(3,j)
   Mtot = Mtot + atomic_weight(atomic_number(j))
 END DO
 
 center_of_mass=sumvector/Mtot
 
 WRITE(output_FID,'(a,3f10.6)')' center_of_mass: ', center_of_mass
 WRITE(output_FID,*)' '
 FLUSH(output_FID)

 END SUBROUTINE compute_center_of_mass

 END MODULE module_compute_center_of_mass
 
 