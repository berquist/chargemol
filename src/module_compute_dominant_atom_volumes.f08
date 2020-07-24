 MODULE module_compute_dominant_atom_volumes
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE
 
 CONTAINS
 
 SUBROUTINE compute_dominant_atom_volumes()
 !===================================================================================
 ! This file assigns each point in the unit cell to the atom with the highest density
 ! (as defined by dominant_atom_weight).
 !===================================================================================
 
 REAL(kind=dp) :: temp_vector(3),distance
 REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: maximum_dominant_atom_weight
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting compute_dominant_atom_volumes'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 ALLOCATE(dominant_atom_points(totnumA,totnumB,totnumC))
 ALLOCATE(maximum_dominant_atom_weight(totnumA,totnumB,totnumC))
 dominant_atom_points=0
 maximum_dominant_atom_weight = 0.0_dp
!$omp parallel default(none) private(na,nb,nc,ka,kb,kc,temp_vector,distance,shell_index,j) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,center_shift,dominant_atom_weight,maximum_dominant_atom_weight,dominant_atom_points,chunk_size,natoms)
 DO j=1,natoms
!$omp do schedule(dynamic,chunk_size)
   DO nc = lower_nc(j),upper_nc(j)
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
           IF (dominant_atom_weight(shell_index,j) > maximum_dominant_atom_weight(ka,kb,kc)) THEN
             dominant_atom_points(ka,kb,kc) = j
           END IF
           maximum_dominant_atom_weight(ka,kb,kc) = max(maximum_dominant_atom_weight(ka,kb,kc),dominant_atom_weight(shell_index,j))
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel 
 !Free some space
 DEALLOCATE(maximum_dominant_atom_weight)
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished compute_dominant_atom_volumes in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE compute_dominant_atom_volumes
 
 END MODULE module_compute_dominant_atom_volumes