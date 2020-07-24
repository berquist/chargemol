 MODULE module_align_periodic_atoms
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
   
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE align_periodic_atoms()
 !===================================================================================
 !Making sure the atomic coordinates for periodic directions are inside the reference
 !unit cell.
 !===================================================================================
  
 INTEGER :: i
  
 WRITE(output_FID,*)'Making sure the atomic coordinates for periodic directions are &
 &inside the reference unit cell.'
 FLUSH(output_FID)
 DO j= 1,natoms
    IF (periodicA) THEN
       center_nabc(1,j) = modulo(center_nabc(1,j),totnumA)
    END IF
    IF (periodicB) THEN
       center_nabc(2,j) = modulo(center_nabc(2,j),totnumB)
    END IF
    IF (periodicC) THEN
       center_nabc(3,j) = modulo(center_nabc(3,j),totnumC)
    END IF
    DO i= 1,3
        coords(i,j) = center_nabc(1,j)*boundary(1,i) + center_nabc(2,j)*boundary(2,i) + &
        center_nabc(3,j)*boundary(3,i) + center_shift(i,j) + origin(i)
    END DO
 END DO
 WRITE(output_FID,*)'The adjusted atomic coordinates are'
 FLUSH(output_FID)
 DO i=1,natoms
   WRITE(output_FID, '(3f10.4)') (coords(j,i),j=1,3)
 END DO
 FLUSH(output_FID)
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'The adjusted center_nabc is'
 FLUSH(output_FID)
 DO i=1,natoms
   WRITE(output_FID, *) (center_nabc(j,i),j=1,3)
 END DO
 FLUSH(output_FID)
  
 END SUBROUTINE align_periodic_atoms 

 END MODULE module_align_periodic_atoms