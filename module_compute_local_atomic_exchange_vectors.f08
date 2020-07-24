 MODULE module_compute_local_atomic_exchange_vectors
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.
 
 USE module_common_variable_declarations
 USE module_global_parameters
 USE module_spin_functions

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE compute_local_atomic_exchange_vectors()
 !=================================================================================== 
 
 INTEGER :: k,comp
 REAL(kind=dp) :: distance,local_density,spin_mag,temp
 REAL(kind=dp),DIMENSION(3) :: temp_vector
 INTEGER,ALLOCATABLE, DIMENSION(:) :: temp_sum_points,temp_sum_density
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting compute_local_atomic_exchange_vectors'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 ALLOCATE(local_spherical_avg_atomic_exchange_vectors(num_nonzero_components,nshells,natoms))
 DO j=1,natoms
   DO k=1,nshells
     local_spherical_avg_atomic_exchange_vectors(1,k,j) = corrected_spherical_avg_density(k,j)
     IF (spin_available) THEN
       IF (non_collinear) THEN
         DO comp=1,3
           local_spherical_avg_atomic_exchange_vectors((comp+1),k,j) = spherical_average_atomic_spin_vector(comp,k,j) 
         END DO  
       ELSE
         local_spherical_avg_atomic_exchange_vectors(2,k,j) = spherical_average_atomic_spin(k,j) 
       END iF    
     END IF
   END DO
 END DO    
 !Sum the local_spherical_avg_atomic_exchange_vectors at each unit cell grid point
 ALLOCATE(total_local_spherical_avg_atomic_exchange_vectors(num_nonzero_components,totnumA,totnumB,totnumC))
 total_local_spherical_avg_atomic_exchange_vectors = 0.0_dp
!$omp parallel default(none) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,center_shift,num_nonzero_components,total_local_spherical_avg_atomic_exchange_vectors,&
!$omp local_spherical_avg_atomic_exchange_vectors,natoms,chunk_size)&
!$omp private(j,nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index,comp)
 DO j=1,natoms
!$omp do &
!$omp schedule(static,chunk_size)
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
           DO comp=1,num_nonzero_components
             total_local_spherical_avg_atomic_exchange_vectors(comp,ka,kb,kc) = total_local_spherical_avg_atomic_exchange_vectors&
             (comp,ka,kb,kc) + local_spherical_avg_atomic_exchange_vectors(comp,shell_index,j)
           END DO
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel

 !Compute the min_density for each atomic radial shell
 ALLOCATE(min_density(nshells,natoms)) 
 min_density = 1000000000.0_dp
!$omp parallel default(none) &
!$omp shared(lower_nc,upper_nc,kCpoints,delta_nc,lower_nb,upper_nb,kBpoints,delta_nb,lower_na,upper_na,kApoints,delta_na,&
!$omp boundary,center_shift,total_local_spherical_avg_atomic_exchange_vectors,natoms,chunk_size)&
!$omp private(j,nc,kc,nb,kb,na,ka,temp_vector,distance,shell_index) &
!$omp reduction(min:min_density)
 DO j=1,natoms
!$omp do &
!$omp schedule(static,chunk_size)
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
           min_density(shell_index,j) = min(min_density(shell_index,j),&
           total_local_spherical_avg_atomic_exchange_vectors(1,ka,kb,kc))
         END IF
       END DO    
     END DO
   END DO
!$omp end do
 END DO
!$omp end parallel

 !Compute the dot product of local total_spherical_avg_atomic_exchange_vectors
 ALLOCATE(dot_product_total_spherical_avg_atomic_exchange_vectors(totnumA,totnumB,totnumC))
 dot_product_total_spherical_avg_atomic_exchange_vectors = 0.0_dp
!$omp parallel do default(none) &
!$omp shared(totnumC,totnumB,totnumA,num_nonzero_components,dot_product_total_spherical_avg_atomic_exchange_vectors,&
!$omp total_local_spherical_avg_atomic_exchange_vectors)&
!$omp private(j,kc,kb,ka,comp)&
!$omp schedule(static,1)
 DO kc = 1,totnumC
   DO kb = 1,totnumB
     DO ka = 1,totnumA
       DO comp=1,num_nonzero_components
         dot_product_total_spherical_avg_atomic_exchange_vectors(ka,kb,kc)=dot_product_total_spherical_avg_atomic_exchange_vectors&
         (ka,kb,kc) + total_local_spherical_avg_atomic_exchange_vectors(comp,ka,kb,kc)**2
       END DO
     END DO
   END DO
 END DO
 !$omp end parallel do
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished compute_local_atomic_exchange_vectors in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 END SUBROUTINE compute_local_atomic_exchange_vectors
 
 END MODULE module_compute_local_atomic_exchange_vectors
