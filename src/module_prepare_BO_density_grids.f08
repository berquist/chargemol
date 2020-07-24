 MODULE module_prepare_BO_density_grids
 ! Copyright (c) 2014, 2015, 2016, 2017 by Thomas A. Manz and Nidia Gabaldon Limas. Rights reserved.

 USE module_common_variable_declarations
 USE module_global_parameters

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE prepare_BO_density_grids()
 !===================================================================================
 
 INTEGER :: k
 REAL(kind=dp) :: distance
 REAL(kind=dp),ALLOCATABLE,DIMENSION(:) :: nA_correction,sum_cubes,temp_sum_density
 REAL(kind=dp),DIMENSION(3) :: temp_vector
 INTEGER,ALLOCATABLE,DIMENSION(:) :: temp_sum_points
 
 INTEGER :: clock_count,clock_count_rate,clock_count_max,clock_count2
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 WRITE(output_FID,*)'Starting prepare_BO_density_grids'
 FLUSH(output_FID)
 CALL system_clock(clock_count,clock_count_rate,clock_count_max)
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 

 WRITE(output_FID,*)'Preparing density grids for bond order analysis ...'
 !Initialize arrays
 ALLOCATE(nA_correction(natoms))
 nA_correction = 0.0_dp
 ALLOCATE(corrected_total_density(totnumA,totnumB,totnumC))
 corrected_total_density = 0.0_dp
 corrected_total_density = total_density
 DEALLOCATE(total_density)
 ALLOCATE(corrected_spherical_avg_density(nshells,natoms))
 corrected_spherical_avg_density = spherical_average_density !Initial estimate
 ALLOCATE(integrated_nA(natoms))
 integrated_nA = 0.0_dp
 ALLOCATE(sum_cubes(natoms))
 ALLOCATE(temp_sum_density(nshells))
 ALLOCATE(temp_sum_points(nshells))
 sum_cubes = 0.0_dp
 K_factor = 0.0_dp
 max_density = 0.0_dp
 !Perform the correction
 DO i = 1,40
   !Compute the correction size
   DO j=1,natoms
     integrated_nA(j) = 0.0_dp
     max_density(j) = 0.0_dp
     DO k=1,nshells
       integrated_nA(j) = integrated_nA(j) + corrected_spherical_avg_density(k,j)*sum_points(k,j)*pixelvolume
     END DO
     nA_correction(j) = valence_population(j) + core_electrons(j) - integrated_nA(j)
     sum_cubes(j) = 0.0_dp
     DO nc = -2,2
       IF (periodicC) THEN
         kc = modulo((nc + center_nabc(3,j)),totnumC) + 1
       ELSE
         kc = nc + center_nabc(3,j) + 1
       END IF
       IF ((kc < 1) .or. (kc > totnumC)) THEN
         CYCLE
       END IF
       DO nb = -2,2
         IF (periodicB) THEN
           kb = modulo((nb + center_nabc(2,j)),totnumB) + 1
         ELSE
           kb = nb + center_nabc(2,j) + 1
         END IF
         IF ((kb < 1) .or. (kb > totnumB)) THEN
           CYCLE
         END IF
         DO na = -2,2
           IF (periodicA) THEN
             ka = modulo((na + center_nabc(1,j)),totnumA) + 1
           ELSE
             ka = na + center_nabc(1,j) + 1
           END IF
           IF ((ka < 1) .or. (ka > totnumA)) THEN
             CYCLE
           END IF
           sum_cubes(j) = sum_cubes(j) + corrected_total_density(ka,kb,kc)**3
           max_density(j) = max(max_density(j),corrected_total_density(ka,kb,kc))
         END DO
       END DO
     END DO
     !Compute the adjustment constant
     K_factor(j) = nA_correction(j)/(sum_cubes(j)*pixelvolume)
     IF (max_density(j) > zero_tolerance) THEN
       K_factor(j) = min(K_factor(j), 0.25_dp/(max_density(j)**2))
     END IF       
     DO nc = -2,2
       IF (periodicC) THEN
         kc = modulo((nc + center_nabc(3,j)),totnumC) + 1
       ELSE
         kc = nc + center_nabc(3,j) + 1
       END IF
       IF ((kc < 1) .or. (kc > totnumC)) THEN
         CYCLE
       END IF
       DO nb = -2,2
         IF (periodicB) THEN
           kb = modulo((nb + center_nabc(2,j)),totnumB) + 1
         ELSE
           kb = nb + center_nabc(2,j) + 1
         END IF
         IF ((kb < 1) .or. (kb > totnumB)) THEN
           CYCLE
         END IF
         DO na = -2,2
           IF (periodicA) THEN
             ka = modulo((na + center_nabc(1,j)),totnumA) + 1
           ELSE
             ka = na + center_nabc(1,j) + 1
           END IF
           IF ((ka < 1) .or. (ka > totnumA)) THEN
             CYCLE
           END IF
           IF (corrected_total_density(ka,kb,kc) > zero_tolerance) THEN
             corrected_total_density(ka,kb,kc) = corrected_total_density(ka,kb,kc)/sqrt(1.0_dp - 2.0_dp*K_factor(j)*&
             corrected_total_density(ka,kb,kc)**2)
           END IF
         END DO
       END DO
     END DO                 
   END DO
   sum_points = 0.0_dp
   sum_density = 0.0_dp
   DO j=1,natoms
     temp_sum_density = 0.0_dp
     temp_sum_points = 0
!$omp parallel do default(none) &
!$omp shared(lower_nc,upper_nc,lower_nb,upper_nb,lower_na,upper_na,j,delta_nc,delta_nb,delta_na,boundary,center_shift,&
!$omp total_pseudodensity,corrected_total_density,partial_density,kCpoints,kBpoints,kApoints,chunk_size) &
!$omp private(nc,nb,na,kc,kb,ka,temp_vector,distance,shell_index) &
!$omp schedule(static,chunk_size) &
!$omp reduction(+:temp_sum_density,temp_sum_points)
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
           IF ((shell_index <= nshells) .and. (total_pseudodensity(ka,kb,kc) > zero_tolerance) ) THEN
             temp_sum_density(shell_index) = temp_sum_density(shell_index) + corrected_total_density(ka,kb,kc)*partial_density(&
             shell_index,j)/total_pseudodensity(ka,kb,kc)
             temp_sum_points(shell_index) = temp_sum_points(shell_index) + 1
           END IF
         END DO    
       END DO
     END DO
!$omp end parallel do
   sum_density(:,j) = temp_sum_density(:)
   sum_points(:,j) = temp_sum_points(:)
   END DO 
   DO j = 1,natoms
     DO k = 1,nshells
       IF (sum_points(k,j) > 0) THEN
         corrected_spherical_avg_density(k,j) = sum_density(k,j)/sum_points(k,j)
       ELSE 
         corrected_spherical_avg_density(k,j) = 0.0_dp
       END IF
     END DO
   END DO
   IF (sum(abs(nA_correction)) < charge_convergence_tolerance) THEN
     WRITE(output_FID,'(A,i5)')' Total density grid corrected in the following number of iterations: ',i
     WRITE(output_FID,'(A,f15.6)')' The total integrated number of electrons is ',sum(sum_density)*pixelvolume
     EXIT
   END IF
 END DO
 
 WRITE(output_FID,*)'Grid preparation complete.'
 
 WRITE(output_FID,*)' '
 WRITE(output_FID,*)'******************** TIME CONTROL ********************'
 CALL system_clock(clock_count2,clock_count_rate,clock_count_max)
 IF (clock_count_rate == 0) THEN
   WRITE(output_FID,*)'Processor clock not available, the time calculation wil not be performed.'
 ELSE
   seconds = modulo((clock_count2-clock_count),clock_count_max)/real(clock_count_rate)
   WRITE(output_FID,*)'Finished prepare_BO_density_grids in ',seconds, 'seconds'
   FLUSH(output_FID)
 END IF
 WRITE(output_FID,*)'******************************************************'
 WRITE(output_FID,*)' '
 FLUSH(output_FID) 
 
 END SUBROUTINE prepare_BO_density_grids

 END MODULE module_prepare_BO_density_grids