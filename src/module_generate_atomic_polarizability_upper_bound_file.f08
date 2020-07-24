 MODULE module_generate_atomic_polarizability_upper_bound_file
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_atomic_number_to_symbol
 USE module_string_utilities

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE generate_atomic_polarizability_upper_bound_file( )
 !===================================================================================

 INTEGER :: polarizability_file,z
 REAL(kind=dp) :: V1(3), V2(3),V3(3),temp_variable
 CHARACTER(45) :: filename
 CHARACTER(2) :: atomic_symbol

 !Print atomic polarizabilities information
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 filename = 'DDEC_atomic_polarizabilities_upper_bound.xyz'
 OPEN(NEWUNIT=polarizability_file, FILE=filename, STATUS='REPLACE')
 WRITE(polarizability_file,'(i5)') natoms
 
 IF (periodicA .or. periodicB .or. periodicC) THEN
   WRITE(polarizability_file,'(a,3f10.6,a,3f10.6,a,3f10.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell [{ ',V1(1),&
   V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'
   FLUSH(polarizability_file)
 ELSE    
    WRITE(polarizability_file,'(a18)')'Nonperiodic system'
 END IF
 
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
   WRITE(polarizability_file,'(a,4f11.6)') atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/bohrperangstrom,&
   final_result(i,5)/bohrperangstrom,atomic_polarizability_upper_bound(i) 
   FLUSH(polarizability_file)
 END DO
 WRITE(polarizability_file,*)' '
  
 ! Print information at the bottom of the file
 WRITE(polarizability_file,'(a)')'Same information as above printed with atom number.'
 DO i = 1,natoms
   z = nint(final_result(i,2))
   CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
   WRITE(polarizability_file,'(i5,a,a,4f11.6)') i,' ',atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/&
   bohrperangstrom,final_result(i,5)/bohrperangstrom,atomic_polarizability_upper_bound(i) 
   FLUSH(polarizability_file)
 END DO
 WRITE(polarizability_file,*)' '
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)
  
 WRITE(polarizability_file,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),":",time(5:6)
 FLUSH(polarizability_file)
 
 CLOSE(polarizability_file)
 
 END SUBROUTINE generate_atomic_polarizability_upper_bound_file
 
 END MODULE module_generate_atomic_polarizability_upper_bound_file