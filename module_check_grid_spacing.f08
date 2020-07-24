 MODULE module_check_grid_spacing
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_precision
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_generate_kApoints_kBpoints_kCpoints

 IMPLICIT NONE
 
 CONTAINS 
 
 SUBROUTINE check_grid_spacing()
 !===================================================================================
 ! Check the minimum pixel volume.
 !===================================================================================
 
 REAL(kind=dp) :: distance,temp_vector(3)
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: population
 INTEGER :: i,j,k
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting check_grid_spacing'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 IF (pixelvolume > maxpixelvolume) THEN
   WRITE(output_FID,*) 'The grid spacing is too coarse.' 
   FLUSH(output_FID)
   WRITE(output_FID,*) 'Please correct your input files and re-submit. Program will terminate.'
   FLUSH(output_FID)
   STOP
 END IF
 !Initialize the partial density array with a guess
 ALLOCATE(partial_density(nshells,natoms))
 partial_density=0.0_dp
 DO i = 1,natoms
   DO j = 1,nshells
        partial_density(j,i) = 0.05_dp*((1.0_dp - (j-1.0_dp)/(nshells - 1.0_dp)))**3
   END DO
 END DO
 CALL generate_kApoints_kBpoints_kCpoints
 !Check to see whether the grid spacing is adequate
 ALLOCATE(population(natoms))
 ALLOCATE(sum_points(nshells,natoms))
 ALLOCATE(sum_density(nshells,natoms))
 population=0.0_dp
 sum_density=0.0_dp
 sum_points = 0
!$omp parallel default(none) private(j,na,nb,nc,temp_vector,distance,shell_index) &
!$omp shared(natoms,lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,center_shift,partial_density,chunk_size) &
!$omp reduction(+:sum_density,sum_points)
 DO j=1,natoms
!$omp do &
!$omp schedule(dynamic,chunk_size)
   DO nc = lower_nc(j),upper_nc(j)
     DO nb = lower_nb(j),upper_nb(j)
       DO na = lower_na(j),upper_na(j)
         temp_vector(1) = na*boundary(1,1) + nb*boundary(2,1) + nc*boundary(3,1) - center_shift(1,j)
         temp_vector(2) = na*boundary(1,2) + nb*boundary(2,2) + nc*boundary(3,2) - center_shift(2,j)
         temp_vector(3) = na*boundary(1,3) + nb*boundary(2,3) + nc*boundary(3,3) - center_shift(3,j)
         distance = sqrt(temp_vector(1)*temp_vector(1) + temp_vector(2)*temp_vector(2) + temp_vector(3)*temp_vector(3))
         shell_index = ceiling(scalefactor*distance + zero_tolerance)
         IF (shell_index <= nshells) THEN
           sum_density(shell_index,j) = sum_density(shell_index,j) + partial_density(shell_index,j)
           sum_points(shell_index,j) = sum_points(shell_index,j) + 1
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO 
!$omp end parallel
 population=sum(sum_density,1)*pixelvolume
 DO j = 1,natoms
   IF (population(j) - sum(population)/natoms > integration_tolerance) THEN
     WRITE(output_FID,*) 'Integration volumes are not sufficiently accurate.'
     FLUSH(output_FID)
     WRITE(output_FID,*) 'Either your electron density input file(s) are too coarsely grained' 
     FLUSH(output_FID)
     WRITE(output_FID,*) 'or they do not include large enough distance(s) along the nonperiodic direction(s) (if any).' 
     FLUSH(output_FID)
     WRITE(output_FID,*) 'Please correct your input files and re-submit. Program will terminate.'
     FLUSH(output_FID)
     WRITE(output_FID,*) 'sum_density'
     FLUSH(output_FID)
     DO i=1,natoms
       WRITE(output_FID,'(10f13.6)') (sum_density(k,i),k=1,nshells)
       FLUSH(output_FID)
     END DO
     WRITE(output_FID,*) 'sum_points'
     FLUSH(output_FID)
     DO i=1,natoms
       WRITE(output_FID,*) (sum_points(k,i),k=1,nshells)
       FLUSH(output_FID)
     END DO
     STOP
   END IF
 END DO
 WRITE(output_FID,*) 'The grid spacing in your electron density input file is adequate.'
 FLUSH(output_FID)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished check_grid_spacing in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE check_grid_spacing
 
 END MODULE module_check_grid_spacing
