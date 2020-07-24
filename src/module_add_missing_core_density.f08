 MODULE module_add_missing_core_density
! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_charge_center_positions_and_parallelpiped
 USE module_check_grid_spacing
 USE module_string_utilities

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE add_missing_core_density( )
 !===================================================================================
  
 REAL(kind=dp), ALLOCATABLE, DIMENSION (:,:) :: ref_core_density
 REAL(kind=dp) :: temp_vector(3),distance
 INTEGER :: k,combined_string_fid,openstat,integer_cutoff
 CHARACTER (200) :: combinedstring
 CHARACTER (200) :: line_of_text
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting add_missing_core_density'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 !First determine the atomic positions and check the grid spacing
 CALL charge_center_positions()
 CALL parallelpiped()
 CALL check_grid_spacing() 
 WRITE(output_FID,*)'Checking to see that all core electrons are accounted for:'
 FLUSH(output_FID)
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'Printing atom number, atomic number, core electrons, and missing core for each atom.'
 WRITE(output_FID,*)'(Missing core electrons will be inserted using stored core electron reference densities.)'
 ALLOCATE(ref_core_density(nshells,natoms))
 integer_cutoff=nint(cutoff_radius)
 DO j = 1,natoms
   WRITE(output_FID,'(2i5,f13.6,i5)')j, atomic_number(j), core_electrons(j), missing_core(j) 
 END DO
 FLUSH(output_FID)
 DO j = 1,natoms
   IF (missing_core(j) == 0) THEN
     CYCLE
   END IF
   !Construct the file name to read
   WRITE(combinedstring, '(A,I3.3,A,I3.3,A,I3.3,A,I3.3,A,I3.3,A)') 'core_',atomic_number(j),'_',atomic_number(j),'_',&
   &missing_core(j),'_',integer_cutoff,'_',nshells,'.txt'
   OPEN(NEWUNIT=combined_string_fid, FILE=(TRIM(ADJUSTL(atomic_densities_directory)))//combinedstring, STATUS='OLD',&
   IOSTAT=openstat)
   IF (openstat > 0) THEN  
     WRITE(output_FID,*) combinedstring
     FLUSH(output_FID)
     WRITE(output_FID,*) "There was a problem reading the atomic density file listed above. Program will terminate."
     FLUSH(output_FID)
     STOP
   END IF
   !Read in the density
   DO
     READ(combined_string_fid,'(a)',IOSTAT=iostat_value) line_of_text
     IF (iostat_value > 0) THEN
       WRITE(output_FID,*) "Could not read the atomic density. Program will terminate."
       FLUSH(output_FID)
       STOP
     ELSE IF (iostat_value < 0) THEN
       EXIT
     ELSE
       IF(line_of_text == 'density in atomic units:')THEN
         DO k=1,nshells
           READ(combined_string_fid,*)ref_core_density(k,j)
         END DO
       END IF
     END IF
   END DO
   CLOSE(combined_string_fid)
 END DO
!$omp parallel default (none) &
!$omp private(j,na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,center_shift,core_density,ref_core_density,chunk_size,natoms,missing_core)
 DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
   DO nc = lower_nc(j),upper_nc(j)
     IF (missing_core(j) > 0) THEN
       kc = kCpoints(j,delta_nc + nc + 1)
       DO nb = lower_nb(j),upper_nb(j)
         kb = kBpoints(j,delta_nb + nb + 1)
         DO na = lower_na(j),upper_na(j)
           ka = kApoints(j,delta_na + na + 1)
           temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
           temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
           temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
           distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
           shell_index = ceiling(scalefactor*distance + zero_tolerance)
           IF (shell_index <= nshells) THEN
             core_density(ka,kb,kc) = core_density(ka,kb,kc) + ref_core_density(shell_index,j)
           END IF
         END DO    
       END DO
     END IF
   END DO
!$omp end do
 END DO
!$omp end parallel
 WRITE(output_FID,*)'Finished the check for missing core electrons.' 
 FLUSH(output_FID)


 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished add_missing_core_density in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE add_missing_core_density
 
 END MODULE module_add_missing_core_density