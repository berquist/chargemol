 MODULE module_generate_atomic_spin_moment_file
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_atomic_number_to_symbol
 USE module_string_utilities 

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE generate_atomic_spin_moment_file( )
 !===================================================================================

 INTEGER :: spin_file_FID,z
 REAL(kind=dp) :: V1(3), V2(3),V3(3),temp_variable
 CHARACTER(50) :: filename
 CHARACTER(2) :: atomic_symbol
 
 V1 = vector1*periA/bohrperangstrom
 V2 = vector2*periB/bohrperangstrom
 V3 = vector3*periC/bohrperangstrom
 IF (StrLowCase(trim(adjustl(charge_type))) == 'ddec3') THEN
   filename = 'DDEC3_atomic_spin_moments.xyz'
 ELSE
   filename = 'DDEC6_even_tempered_atomic_spin_moments.xyz'
 END IF 
 OPEN(NEWUNIT=spin_file_FID, FILE=filename, STATUS='REPLACE')
 WRITE(spin_file_FID,'(i5)') natoms
 
 IF (periodicA .or. periodicB .or. periodicC) THEN
   WRITE(spin_file_FID,'(a,3f10.6,a,3f10.6,a,3f10.6,a)') 'jmolscript: load "" {1 1 1} spacegroup "x,y,z" unitcell [{ ',V1(1),&
   V1(2),V1(3),' }, { ',V2(1),V2(2),V2(3),' }, { ',V3(1),V3(2),V3(3),' }]'
   FLUSH(spin_file_FID)
 ELSE    
    WRITE(spin_file_FID,'(a18)')'Nonperiodic system'
 END IF
 
 DO i = 1,natoms
   IF (spin_available) THEN
     IF (non_collinear) THEN  
       z = nint(final_result(i,2))
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       WRITE(spin_file_FID,'(a,7f11.6)') atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/bohrperangstrom,&
       final_result(i,5)/bohrperangstrom,sqrt(spin_population_vector(1,i)**2 + spin_population_vector(2,i)**2 + &
       &spin_population_vector(3,i)**2),spin_population_vector(1,i),spin_population_vector(2,i),spin_population_vector(3,i)
       FLUSH(output_FID)
     ELSE
       z = nint(final_result(i,2))
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       WRITE(spin_file_FID,'(a,4f11.6)') atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/bohrperangstrom,&
       final_result(i,5)/bohrperangstrom,spin_population(i) 
       FLUSH(spin_file_FID)
      END IF
   ELSE
     z = nint(final_result(i,2))
     CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
     WRITE(spin_file_FID,'(a,3f11.6)') atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/bohrperangstrom,&
     final_result(i,5)/bohrperangstrom
     FLUSH(output_FID)
   END IF
 END DO
 ! Print information at the bottom of the file
 WRITE(spin_file_FID,*)'   '  
 FLUSH(spin_file_FID)
 IF (spin_available) THEN
   IF (non_collinear) THEN
     WRITE(spin_file_FID,*)'Noncollinear spin population analysis was performed '  
     FLUSH(spin_file_FID)
     WRITE(spin_file_FID,'(a,3f10.6)')'The total spin magnetic moment vector of the unit cell is  ',tot_spin_moment_vector(1)&
     ,tot_spin_moment_vector(2),tot_spin_moment_vector(3)
     FLUSH(spin_file_FID)
   ELSE
     WRITE(spin_file_FID,'(a48)') 'Collinear spin population analysis was performed '
     FLUSH(spin_file_FID)
     WRITE(spin_file_FID,'(a,f10.6)') 'The total spin magnetic moment of the unit cell is  ',tot_spin_moment
     FLUSH(spin_file_FID)
     WRITE(spin_file_FID,*)'   '   
     FLUSH(spin_file_FID)
   END IF
 END IF

 ! Print information at the bottom of the file
 WRITE(spin_file_FID,*)version   
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a44)')'See ddec.sourceforge.net for latest version.'   
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,*)'   '   
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a25)')'Computational parameters:'   
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a,f4.2)') 'Spin_ref_fraction = ',spin_ref_fraction
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a,i3)') 'Number of radial integration shells = ',nshells
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a,i4)') 'Cutoff radius (pm) = ',nint(cutoff_radius)   
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a,f8.6)') 'Spin convergence tolerance = ',spin_convergence_tolerance 
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,'(a,i2)') 'Number of iterations to convergence = ',spin_iter  
 FLUSH(spin_file_FID)
 WRITE(spin_file_FID,*)'   ' 
 FLUSH(spin_file_FID)
 
 ! Print ASMs with atom number
 DO i = 1,natoms
   IF (.not. spin_available) THEN
     z = nint(final_result(i,2))
     CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
     WRITE(spin_file_FID,'(i5,a,a,3f11.6)') i,' ',atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/&
     bohrperangstrom,final_result(i,5)/bohrperangstrom
     FLUSH(output_FID)
   END IF
   IF (spin_available) THEN
     IF (non_collinear) THEN  
       z = nint(final_result(i,2))
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       WRITE(spin_file_FID,'(i5,a,a,7f11.6)')i,' ',atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/&
       bohrperangstrom,final_result(i,5)/bohrperangstrom,sqrt(spin_population_vector(1,i)**2 + spin_population_vector(2,i)**2 + &
       &spin_population_vector(3,i)**2),spin_population_vector(1,i),spin_population_vector(2,i),spin_population_vector(3,i)
       FLUSH(output_FID)
     ELSE
       z = nint(final_result(i,2))
       CALL atomic_number_to_symbol(z,output_FID,atomic_symbol)
       WRITE(spin_file_FID,'(i5,a,a,4f11.6)') i,' ',atomic_symbol,final_result(i,3)/bohrperangstrom,final_result(i,4)/&
       bohrperangstrom,final_result(i,5)/bohrperangstrom,spin_population(i) 
       FLUSH(spin_file_FID)
      END IF
   END IF
 END DO
 
 CALL date_and_time(date,time,zone,values)
 CALL date_and_time(DATE=date,ZONE=zone)
 CALL date_and_time(TIME=time)
 CALL date_and_time(VALUES=values)

 WRITE(spin_file_FID,*) date(1:4),"/",date(5:6),"/",date(7:),"  ",time(1:2),":",time(3:4),":",time(5:6)
 FLUSH(spin_file_FID)
 
 CLOSE(spin_file_FID)
 
 END SUBROUTINE generate_atomic_spin_moment_file
 
 END MODULE module_generate_atomic_spin_moment_file